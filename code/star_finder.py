import pandas as pd
from pathlib import Path
from pathlib import PurePath
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
image_file = get_pkg_data_filename('20240307_data_I.fits')
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
mean, median, std = sigma_clipped_stats(hdu.data, sigma=3.0)

#hiermee kan je een deel van de data nemen ipv de hele afbeelding
section1 = hdu.data[0:3900,35:4500] #hele afbeelding met ruisgedeelten eruit
plt.subplot(111,projection=wcs) 
mask = np.zeros(section1.shape, dtype=bool)
mask[1620:1850,2425:2630] = True
mask[1800:2050,2730:2975] = True
mask[2180:2320,3920:4080] = True
mask[2020:2070,2685:2700] = True #ruis
mask[965:990,1780:1808] = True
plt.imshow(section1, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 
tbl = find_peaks(section1, median+5*std, box_size=11, mask=mask, border_width=10, wcs=wcs)
tbl['peak_value'].info.format = '%.8g'

tbl['id'] = np.arange(len(tbl))
positions = np.transpose((tbl['x_peak'], tbl['y_peak']))
apertures = CircularAperture(positions, r=25)
apertures.plot(color='red', lw=2, alpha=0.5)

#Dit deel voegt de datapunten toe van sterren die uitgemask waren in het eerste gedeelte
#plt.subplot(222,projection=wcs) 
#plt.imshow(section1, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 
mask1_2 = np.zeros(section1.shape, dtype=bool)
mask1_2[35:4500,0:2000] = True
tbl1_2 = find_peaks(section1, median+5*std, box_size=50, border_width=10, wcs=wcs, npeaks=3, mask=mask1_2)
positions1_2 = np.transpose((tbl1_2['x_peak'], tbl1_2['y_peak']))
apertures1_2 = CircularAperture(positions1_2, r=25)
apertures1_2.plot(color='red', lw=2, alpha=0.5)

#data formatten
tbl1_2['id'] = np.arange(len(tbl1_2))+len(tbl)
for id in range(len(tbl1_2)):
    tbl.add_row(tbl1_2[id])
print(tbl)
plt.show()

#tbl.write('stars_full_img.ecsv', overwrite=True) #file maken met data
