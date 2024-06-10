from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAperture
from photutils.detection import find_peaks

#zet de locatie van member data hier (het zit in hetzelfde mapje dus hoeft niet verandert te worden):
image_file = get_pkg_data_filename('20240306_data_I.fits')
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
mean, median, std = sigma_clipped_stats(hdu.data, sigma=3.0)

#hiermee kan je een deel van de data nemen ipv de hele afbeelding
section1 = hdu.data[0:3900,35:4500] #hele afbeelding met ruisgedeelten eruit
plt.subplot(211) 
mask = np.zeros(section1.shape, dtype=bool)
mask[1600:1920,2350:2660] = True
plt.imshow(section1, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 
tbl = find_peaks(section1, median+5*std, box_size=11, mask=mask, border_width=10, wcs=wcs)
tbl['peak_value'].info.format = '%.8g'

tbl['id'] = np.arange(len(tbl))
positions = np.transpose((tbl['x_peak'], tbl['y_peak']))
apertures = CircularAperture(positions, r=25)
apertures.plot(color='red', lw=2, alpha=0.5)

#Dit deel voegt de datapunten toe van sterren die uitgemask waren in het eerste gedeelte
plt.subplot(212) 
plt.imshow(section1, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 
mask1_2 = np.zeros(section1.shape, dtype=bool)
tbl1_2 = find_peaks(section1, median+5*std, box_size=50, border_width=10, wcs=wcs, npeaks=4, mask=mask1_2)
positions1_2 = np.transpose((tbl1_2['x_peak'], tbl1_2['y_peak']))
apertures1_2 = CircularAperture(positions1_2, r=25)
apertures1_2.plot(color='red', lw=2, alpha=0.5)

#data formatten
tbl1_2['id'] = np.arange(len(tbl1_2))+len(tbl)
for id in range(len(tbl1_2)):
    tbl.add_row(tbl1_2[id])
print(tbl)
plt.show()

#tbl.write('stars_full_img_0306_I.ecsv', overwrite=True) #file maken met data
