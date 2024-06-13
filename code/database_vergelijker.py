import pandas as pd
from pathlib import Path
from pathlib import PurePath
from IPython.display import display
import matplotlib.pyplot as plt
import numpy as np
from photutils.aperture import CircularAperture

#flexible file location
main_loc = Path(__file__).resolve().parent.parent
loc = PurePath(main_loc,"data\\asu.tsv")

#data wordt geimporteerd als pandas Dataframe (https://www.geeksforgeeks.org/python-pandas-dataframe/)
data_db = pd.read_csv(loc, sep=';', comment='#')
data_db = data_db.drop([0,1])
data_db['DE_ICRS'] = data_db['DE_ICRS'].str.replace('+', '')
data_db['DE_ICRS'] = data_db['DE_ICRS'].astype(float)
data_db['RA_ICRS'] = data_db['RA_ICRS'].astype(float)
data_db = data_db.rename(columns={'Plx':'plx'})

#import hier database met meer sterren van Mitchel, dan kunnen we data dubbelchecken en vergelijken
loc_gaudain = PurePath(main_loc,"data\\cantat_gaudain2018.tsv")
data_goudain = pd.read_csv(loc_gaudain, sep='\t', comment='#')
data_goudain = data_goudain.drop([0,1])
data_goudain['DE_ICRS'] = data_goudain['DE_ICRS'].str.replace('+', '')
data_goudain['DE_ICRS'] = data_goudain['DE_ICRS'].astype(float)
data_goudain['RA_ICRS'] = data_goudain['RA_ICRS'].astype(float)

#import hier de twee sterdatabases die we zelf hebben!
loc_stars_06 = PurePath(main_loc,"code\\stars_data_20240306.ecsv")
loc_stars_07 = PurePath(main_loc,"code\\stars_data_20240307.ecsv")
data_stars_06 = pd.read_csv(loc_stars_06, sep=' ', comment='#')
data_stars_07 = pd.read_csv(loc_stars_07, sep=' ', comment='#')
dec_corr_07 = -0.004168
data_stars_07['skycoord_peak.dec'] = data_stars_07['skycoord_peak.dec'] + dec_corr_07

plt.figure(1)
plt.scatter(data_stars_06['skycoord_peak.ra'],data_stars_06['skycoord_peak.dec'],color='red')
plt.scatter(data_stars_07['skycoord_peak.ra'],data_stars_07['skycoord_peak.dec'],color='blue')
plt.scatter(x=data_db['RA_ICRS'], y=data_db['DE_ICRS'], color='purple', alpha=0.8)
plt.scatter(x=data_goudain['RA_ICRS'], y=data_goudain['DE_ICRS'], color='green', alpha=0.8)

positions_06  = data_stars_06[['skycoord_peak.ra', 'skycoord_peak.dec']]
apertures_06 = CircularAperture(positions_06, r=0.0001197)
apertures_06.plot(color='red', lw=2, alpha=1)
positions_07  = data_stars_07[['skycoord_peak.ra', 'skycoord_peak.dec']]
apertures_07 = CircularAperture(positions_07, r=0.0001197)
apertures_07.plot(color='blue', lw=2, alpha=1)
positions_alf  = data_db[['RA_ICRS', 'DE_ICRS']]
apertures_alf = CircularAperture(positions_alf, r=1.11*10**(-4))
apertures_alf.plot(color='purple', lw=2, alpha=1)
positions_goud  = data_goudain[['RA_ICRS', 'DE_ICRS']]
apertures_goud = CircularAperture(positions_goud, r=1.11*10**(-4))
apertures_goud.plot(color='green', lw=2, alpha=1)

#for v in data_stars_06[['skycoord_peak.ra', 'skycoord_peak.dec']]:
ra_search_range = 0.000660*2
dec_search_range = 0.0001907*2
""" for id, row in data_stars_06.iterrows():
    for id_db, row_db in data_db.iterrows():
        if abs(row['skycoord_peak.ra'] - row_db['RA_ICRS']) < ra_search_range and abs(row['skycoord_peak.dec'] - row_db['DE_ICRS']) < dec_search_range:
            new_row = [row['skycoord_peak.ra'], row['skycoord_peak.dec']]
            members_list.append(new_row)

members_array = np.array(members_list)
member_stars = pd.DataFrame(data=members_array, columns=['RA', 'DEC']) """

#dit deel is gegenereerd door copilot gebaseerd op de methodes die hierboven gecomment zijn. Verbeterde efficiency
def db_verg(df, df_db):
    df_db['plx'] = df_db['plx'].astype(float)

    ra_stars = df['skycoord_peak.ra'].values
    dec_stars = df['skycoord_peak.dec'].values
    ra_db = df_db['RA_ICRS'].values
    dec_db = df_db['DE_ICRS'].values
    plx = df_db['plx'].values

    x_peak = df['x_peak'].values
    y_peak = df['y_peak'].values
    peak_value = df['peak_value'].values
    id = df['id'].values
    img = df['img'].values
    

    ra_diff = np.abs(ra_stars[:, None] - ra_db)
    dec_diff = np.abs(dec_stars[:, None] - dec_db)

    matches = (ra_diff < ra_search_range) & (dec_diff < dec_search_range)

    # Extract the matching coordinates
    matched_indices = np.any(matches, axis=1)
    matched_ra = ra_stars[matched_indices]
    matched_dec = dec_stars[matched_indices]
    matched_x = x_peak[matched_indices]
    matched_y = y_peak[matched_indices]
    matched_peak_val = peak_value[matched_indices]
    matched_id = id[matched_indices]
    matched_img = img[matched_indices]
    
    matched_plx = []
    for i, star_match in enumerate(matches):
        if np.any(star_match):
            db_index = np.where(star_match)[0][0]  # Take the first matching index
            matched_plx.append(plx[db_index])

    # Create the resulting DataFrame
    members_array = np.column_stack((matched_ra, matched_dec, matched_x, matched_y, matched_peak_val, matched_id, matched_plx, matched_img))
    return pd.DataFrame(data=members_array, columns=['skycoord_peak.ra', 'skycoord_peak.dec', 'x_peak', 'y_peak', 'peak_value', 'id', 'plx', 'img'])

member_stars_06_alf = db_verg(data_stars_06, data_db)
member_stars_07_alf = db_verg(data_stars_07, data_db)
member_stars_06_goud = db_verg(data_stars_06, data_goudain)
member_stars_07_goud = db_verg(data_stars_07, data_goudain)
member_stars = pd.concat([member_stars_06_alf, member_stars_07_alf, member_stars_06_goud, member_stars_07_goud], ignore_index=True)
positions_memb = member_stars[['skycoord_peak.ra', 'skycoord_peak.dec']]
apertures_memb= CircularAperture(positions_memb, r=0.005)
apertures_memb.plot(color='black', lw=2, alpha=1)

print(member_stars)

#afbeeldingen plotten met member stars
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats

plt.figure(2)
image_file = get_pkg_data_filename('20240306_data_I.fits')
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
mean, median, std = sigma_clipped_stats(hdu.data, sigma=3.0)
plt.subplot(211,projection=wcs) 
plt.imshow(hdu.data, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 

pixel_coords = wcs.wcs_world2pix(member_stars[['skycoord_peak.ra', 'skycoord_peak.dec']].values, 1)

member_stars['x_pixel'] = pixel_coords[:, 0]
member_stars['y_pixel'] = pixel_coords[:, 1]
plt.scatter(member_stars['x_pixel'], member_stars['y_pixel'], facecolors='none', edgecolors='r', s=5)

image_file2 = get_pkg_data_filename('20240307_data_I.fits')
hdu2 = fits.open(image_file2)[0]
wcs2 = WCS(hdu2.header)
mean2, median2, std2 = sigma_clipped_stats(hdu2.data, sigma=3.0)
plt.subplot(212,projection=wcs2) 
plt.imshow(hdu2.data, origin='lower', cmap='Greys', vmin=median2, vmax=median2+5*std2, interpolation='nearest') 

pixel_coords = wcs2.wcs_world2pix(member_stars[['skycoord_peak.ra', 'skycoord_peak.dec']].values, 1)
member_stars['x_pixel'] = pixel_coords[:, 0]
member_stars['y_pixel'] = pixel_coords[:, 1]
plt.scatter(member_stars['x_pixel'], member_stars['y_pixel'], facecolors='none', edgecolors='r', s=5)

#member_stars.to_csv('member_stars.csv') #file maken met data
plt.show()
