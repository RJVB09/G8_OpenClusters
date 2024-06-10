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

section1 = hdu.data
plt.subplot(211) 
mask = np.zeros(section1.shape, dtype=bool)
mask[0:300,0:4500] = True
mask[0:4500,0:200] = True
mask[2840:3030,1425:1625] = True
mask[3215:3425,2600:2820] = True
mask[3207:3225,2650:2675] = False
mask[3852:3864,670:685] = True
mask[2365:2380,495:510] = True
mask[1982:2090,3005:3108] = True
mask[2080:2100,3095:3120] = True
mask[1000:1437,2690:3085] = True
mask[1420:1450,2875:2905] = True
mask[1135:1175,2830:2875] = False
mask[1170:1195,2730:2800] = False
mask[1265:1282,2858:2872] = False
mask[1005:1025,2725:2740] = False
mask[1050:1059,2945:2955] = False
mask[1200:1215,3066:3078] = False
mask[875:1100,800:1060] = True
mask[2111:2115,899:904] = True
mask[2217:2450,2830:3085] = True
mask[2240:2272,2965:3005] = False
mask[2080:2365,1680:2050] = True
mask[2292:2310,1945:1965] = False
mask[240:300,3500:3600] = False
mask[220:265,3480:3550] = True
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
mask1_2[0:4500,0:4500] = True
mask1_2[868:3410,825:3130] = False
tbl1_2 = find_peaks(section1, median+5*std, box_size=50, border_width=10, wcs=wcs, npeaks=10, mask=mask1_2)
positions1_2 = np.transpose((tbl1_2['x_peak'], tbl1_2['y_peak']))
apertures1_2 = CircularAperture(positions1_2, r=25)
apertures1_2.plot(color='red', lw=2, alpha=0.5)

#data formatten
tbl1_2['id'] = np.arange(len(tbl1_2))+len(tbl)
for id in range(len(tbl1_2)):
    tbl.add_row(tbl1_2[id])
print(tbl)
plt.show()

tbl.write('stars_full_img_0306_I.ecsv', overwrite=True) #file maken met data
