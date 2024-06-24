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
loc_stars_06 = PurePath(main_loc,"code\\stars_data_20240306.csv")
loc_stars_07 = PurePath(main_loc,"code\\stars_data_20240307.csv")
data_stars_06 = pd.read_csv(loc_stars_06, sep=',', comment='#')
data_stars_07 = pd.read_csv(loc_stars_07, sep=',', comment='#')
dec_corr_07 = -0.004168
data_stars_07['skycoord_peak.dec'] = data_stars_07['skycoord_peak.dec'] + dec_corr_07

plt.figure(1)
plt.subplot(121)
plt.scatter(data_stars_06['skycoord_peak.ra'],data_stars_06['skycoord_peak.dec'],color='#b2df8a',marker='*')
plt.scatter(data_stars_07['skycoord_peak.ra'],data_stars_07['skycoord_peak.dec'],color='#33a02c',marker='*')
plt.scatter(x=data_db['RA_ICRS'], y=data_db['DE_ICRS'], color='#a6cee3', alpha=0.7,marker='8')
plt.scatter(x=data_goudain['RA_ICRS'], y=data_goudain['DE_ICRS'], color='#1f78b4', alpha=0.7,marker='8')
plt.xlabel('RA (degrees)')
plt.ylabel('DEC (degrees)')
plt.gca().set_aspect('equal')
positions_06  = data_stars_06[['skycoord_peak.ra', 'skycoord_peak.dec']]
apertures_06 = CircularAperture(positions_06, r=0.0001197)
apertures_06.plot(color='#b2df8a', lw=2, alpha=1)
positions_07  = data_stars_07[['skycoord_peak.ra', 'skycoord_peak.dec']]
apertures_07 = CircularAperture(positions_07, r=0.0001197)
apertures_07.plot(color='#33a02c', lw=2, alpha=1)
positions_alf  = data_db[['RA_ICRS', 'DE_ICRS']]
apertures_alf = CircularAperture(positions_alf, r=1.11*10**(-4))
apertures_alf.plot(color='#a6cee3', lw=2, alpha=1)
positions_goud  = data_goudain[['RA_ICRS', 'DE_ICRS']]
apertures_goud = CircularAperture(positions_goud, r=1.11*10**(-4))
apertures_goud.plot(color='#1f78b4', lw=2, alpha=1)
plt.legend(["stars img 06/03", "stars img 07/03", "members Alfonso", "members Goudain"], loc="lower right")
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
    ob = df['ob'].values
    

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
    matched_ob = ob[matched_indices]
    
    matched_plx = []
    for i, star_match in enumerate(matches):
        if np.any(star_match):
            db_index = np.where(star_match)[0][0]  # Take the first matching index
            matched_plx.append(plx[db_index])

    # Create the resulting DataFrame
    members_array = np.column_stack((matched_ra, matched_dec, matched_x, matched_y, matched_peak_val, matched_id, matched_plx, matched_img, matched_ob))
    return pd.DataFrame(data=members_array, columns=['skycoord_peak.ra', 'skycoord_peak.dec', 'x_peak', 'y_peak', 'peak_value', 'id', 'plx', 'img', 'ob'])

member_stars_06_alf = db_verg(data_stars_06, data_db)
member_stars_07_alf = db_verg(data_stars_07, data_db)
member_stars_06_goud = db_verg(data_stars_06, data_goudain)
member_stars_07_goud = db_verg(data_stars_07, data_goudain)
member_stars = pd.concat([member_stars_06_alf, member_stars_07_alf, member_stars_06_goud, member_stars_07_goud], ignore_index=True)
member_stars = member_stars.drop_duplicates(subset=['skycoord_peak.ra', 'skycoord_peak.dec'], keep='first', ignore_index=True) #haalt duplicates er uit
positions_memb = member_stars[['skycoord_peak.dec','skycoord_peak.ra']]
apertures_memb= CircularAperture(positions_memb, r=0.005)
apertures_memb.plot(color='black', lw=2, alpha=1)
plt.rc('axes', titlesize=18)
plt.rc('axes', labelsize=18)
plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)
plt.rc('legend', fontsize=13)
plt.rc('figure', titlesize=18)

#meer duplicates
#duplicates = member_stars[(member_stars['id'] == 516) | (member_stars['id'] == 493) | (member_stars['id'] == 351) | (member_stars['id'] == 317) | (member_stars['id'] == 112) | (member_stars['id'] == 153) | (member_stars['id'] == 9) | (member_stars['id'] == 27)]
#member_stars = pd.concat([member_stars, duplicates[duplicates['img']==7], duplicates[duplicates['img']==7]]).drop_duplicates(keep=False)
duplicates2 = member_stars[(member_stars['id'] == 516) | (member_stars['id'] == 493) | (member_stars['id'] == 351) | (member_stars['id'] == 317) | (member_stars['id'] == 112) | (member_stars['id'] == 153) | (member_stars['id'] == 9) | (member_stars['id'] == 27)]

print(member_stars)
plt.figure(4)
plt.subplot(111)
plt.scatter(data_stars_06['skycoord_peak.dec'],data_stars_06['skycoord_peak.ra'],color='#33a02c',marker='*', s=60, zorder=10)
plt.scatter(y=data_db['RA_ICRS'], x=data_db['DE_ICRS'], color='#a6cee3', alpha=0.7,marker='8', zorder=5)
plt.scatter(y=data_goudain['RA_ICRS'], x=data_goudain['DE_ICRS'], color='#1f78b4', alpha=0.7,marker='8', zorder=5)
plt.gca().set_aspect('equal')
plt.xlim(19.3,20.4)
plt.ylim(129.7,130.4)

apertures_memb= CircularAperture(positions_memb, r=0.015)
apertures_memb.plot(color='red', lw=2, alpha=1, zorder=15)

plt.scatter(data_stars_07['skycoord_peak.dec'],data_stars_07['skycoord_peak.ra'],color='#33a02c',marker='*', zorder=10, s=60)
plt.ylabel('RA (degrees)', fontsize=20)
plt.xlabel('DEC (degrees)', fontsize=20)
plt.legend(["Stars from data","Members Alfonso 2023","Members Cantat-Gaudin 2018","Confirmed members",],loc="lower left")


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
plt.subplot(122,projection=wcs) 
plt.imshow(hdu.data, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 
pixel_coords = wcs.wcs_world2pix(member_stars[['skycoord_peak.ra', 'skycoord_peak.dec']].values, 1)

member_stars['x_pixel'] = pixel_coords[:, 0]
member_stars['y_pixel'] = pixel_coords[:, 1]
plt.scatter(member_stars['x_pixel'], member_stars['y_pixel'], facecolors='none', edgecolors='r', s=30)
plt.tick_params(axis='x', which='both', labelbottom=True, labeltop=False)
plt.tick_params(axis='y', which='both', labelbottom=True, labeltop=False)
plt.xlim(0,4100)

image_file2 = get_pkg_data_filename('20240307_data_I.fits')
hdu2 = fits.open(image_file2)[0]
wcs2 = WCS(hdu2.header)
mean2, median2, std2 = sigma_clipped_stats(hdu2.data, sigma=3.0)
plt.subplot(121,projection=wcs2) 
plt.imshow(hdu2.data, origin='lower', cmap='Greys', vmin=median2, vmax=median2+5*std2, interpolation='nearest') 

pixel_coords = wcs2.wcs_world2pix(member_stars[['skycoord_peak.ra', 'skycoord_peak.dec']].values, 1)
member_stars['x_pixel'] = pixel_coords[:, 0]
member_stars['y_pixel'] = pixel_coords[:, 1]
plt.scatter(member_stars['x_pixel'], member_stars['y_pixel'], facecolors='none', edgecolors='r', s=30)
plt.tick_params(axis='x', which='both', labelbottom=True, labeltop=False)
plt.tick_params(axis='y', which='both', labelbottom=True, labeltop=False)
plt.xlim(0,4100)
for index, row in member_stars.iterrows():
    plt.text(row['x_pixel'], row['y_pixel'], int(row['id']), color='r', fontsize='medium')

plt.figure(3)
plt.scatter(member_stars['x_pixel'], member_stars['y_pixel'],facecolors='red', s=5)
for index, row in member_stars.iterrows():
    plt.text(row['x_pixel'], row['y_pixel'], int(row['id']), color='r', fontsize='medium')

#member_stars.to_csv('member_stars_test.csv') #file maken met data

plt.show()
