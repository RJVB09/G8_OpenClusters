import pandas as pd
from pathlib import Path
from pathlib import PurePath

#flexible file location
loc = Path(__file__).resolve().parent.parent
loc = PurePath(loc,"data\\asu.tsv")
print(loc)

#data wordt geimporteerd als pandas Dataframe (https://www.geeksforgeeks.org/python-pandas-dataframe/)
data = pd.read_csv(loc, sep=';', comment='#')
print(data)

#als we een bepaalde ster zoeken uit onze data
RA_low = 0.0
RA_high = 360.0
DE_low = 0.0
DE_high = 360.0

#elk datapunt wordt geprint
for index, row in data.iterrows():
    #index, gaia_id, Right Asc, Declination, Parallax, proper motion RA, proper motion 
    print(index, row['GaiaDR3'], row['RA_ICRS'], row['DE_ICRS'], row['Plx'], row['pmRA'], row['pmDE'])
    #if index > 1 and float(row['RA_ICRS']) > RA_low and float(row['RA_ICRS']) < RA_high and float(row['DE_ICRS']) > DE_low and float(row['DE_ICRS']) < DE_high:
        #print('BINGO') 


#deel 2 :)
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
from matplotlib.colors import LogNorm
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture
from photutils.detection import find_peaks

#zet de locatie van member data hier (het zit in hetzelfde mapje dus hoeft niet verandert te worden):
image_file = get_pkg_data_filename('20240306_data.fits')
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
mean, median, std = sigma_clipped_stats(hdu.data, sigma=3.0)
print('full img stats', mean, median, std)

#hiermee kan je een deel van de data nemen ipv de hele afbeelding
section1 = hdu.data[0:3900,35:4500]
section2 = hdu.data[1700:2400,2250:4500]

#peak techniek
#projection=wcs
plt.subplot(121) 
mask = np.zeros(section1.shape, dtype=bool)
mask[1620:1850,2425:2630] = True
mask[1800:2050,2730:2975] = True
mask[2180:2320,3920:4080] = True
mask[2020:2070,2685:2700] = True #ruis
mask[965:990,1780:1808] = True
plt.imshow(section1, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 
tbl = find_peaks(section1, median+5*std, box_size=11, mask=mask, border_width=10, wcs=wcs)

tbl['peak_value'].info.format = '%.8g'


print(len(tbl))  # print only the first 10 peaks
tbl['id'] = np.arange(len(tbl))
positions = np.transpose((tbl['x_peak'], tbl['y_peak']))
apertures = CircularAperture(positions, r=25)
apertures.plot(color='red', lw=2, alpha=0.5)
print(tbl)

#Dit deel voegt de datapunten toe van sterren die uitgemask waren in het eerste gedeelte
plt.subplot(122) 
plt.imshow(section1, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 
mask2 = np.zeros(section1.shape, dtype=bool)
mask2[35:4500,0:2000] = True
tbl2 = find_peaks(section1, median+5*std, box_size=50, border_width=10, wcs=wcs, npeaks=3, mask=mask2)
positions2 = np.transpose((tbl2['x_peak'], tbl2['y_peak']))

apertures2 = CircularAperture(positions2, r=25)
apertures2.plot(color='red', lw=2, alpha=0.5)

#data formatten
tbl2['id'] = np.arange(len(tbl2))+len(tbl)
for id in range(len(tbl2)):
    tbl.add_row(tbl2[id])


#tbl.write('stars_full_img.ecsv', overwrite=True) #file maken met data
plt.show()

#deel 3 :) DIT IS VOOR DE ANDERE AFBEELDING!