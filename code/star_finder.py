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
image_file = get_pkg_data_filename('test.fits')
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
mean, median, std = sigma_clipped_stats(hdu.data, sigma=3.0)
print('full img stats', mean, median, std)

#hiermee kan je een deel van de data nemen ipv de hele afbeelding
section1 = hdu.data[0:3900,35:4500]


plt.subplot(131, projection=wcs)  
plt.imshow(section1, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 

plt.subplot(132, projection=wcs) 
plt.imshow(section1, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 
#plt.colorbar()

#WORK IN PROCESS :)

#DAO techniek
daofind = DAOStarFinder(fwhm=4.0, threshold=5*std)
sources = daofind(section1-median)
for col in sources.colnames:
    if col not in ('id','npix'):
        sources[col].info.format = '%.2f'
sources.pprint(max_width=200)
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=25)
apertures.plot(color='red', lw=2, alpha=0.5)

#peak techniek
plt.subplot(133, projection=wcs) 
plt.imshow(section1, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 

tbl = find_peaks(section1, median+5*std, box_size=11)
tbl['peak_value'].info.format = '%.8g'  # for consistent table output
print(len(tbl))  # print only the first 10 peaks
positions2 = np.transpose((tbl['x_peak'], tbl['y_peak']))
apertures2 = CircularAperture(positions2, r=25)
apertures2.plot(color='red', lw=2, alpha=0.5)


plt.show()

