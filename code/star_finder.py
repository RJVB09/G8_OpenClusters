import pandas as pd

#zet de locatie van member data hier:
loc = "C:\\python programmas\\G8_OpenClusters\\data\\asu.tsv"

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

#zet de locatie van member data hier:
image_file = get_pkg_data_filename('test.fits')
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
mean, median, std = sigma_clipped_stats(hdu.data, sigma=3.0)
print('full img stats', mean, median, std)

""" plt.subplot(133, projection=wcs)
section1 = hdu.data[0:4500, 3000:4500]
plt.imshow(section1, origin='lower', cmap='Greys', vmin=mean, vmax=mean+5*std)  """


#plt.subplot(131, projection=wcs) 
plt.subplot(projection=wcs) 
plt.imshow(hdu.data, origin='lower', cmap='Greys', vmin=np.mean(hdu.data), vmax=0.0005*np.max(hdu.data)+np.mean(hdu.data)) 

""" plt.subplot(132, projection=wcs) 
plt.imshow(hdu.data, origin='lower', cmap='Greys', vmin=mean, vmax=mean+5*std) 
#plt.colorbar() """

#WORK IN PROCESS :)

daofind = DAOStarFinder(fwhm=3.0, threshold=5*std)
sources = daofind(hdu.data-median)
for col in sources.colnames:
    if col not in ('id','npix'):
        sources[col].info.format = '%.2f'
sources.pprint(max_width=76)
plt.show()

