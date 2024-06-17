from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAperture
from photutils.detection import find_peaks

#file van afbeelding wordt geopend en WCS data extract
image_file = get_pkg_data_filename('20240307_data_I.fits')
hdu = fits.open(image_file)[0]
wcs = WCS(hdu.header)
mean, median, std = sigma_clipped_stats(hdu.data, sigma=3.0)

#heel veel masking haha help...
section1 = hdu.data[0:3900,35:4500] #hele afbeelding met ruisgedeelten eruit
plt.subplot(211) 
mask = np.zeros(section1.shape, dtype=bool)
mask[1600:1920,2350:2660] = True
mask[1750:1920,2620:2760] = True
mask[1780:2050,2730:2975] = True
mask[2180:2320,3920:4080] = True
mask[2020:2070,2685:2700] = True #ruis
mask[965:990,1780:1808] = True
mask[0:25,0:4500] = True
mask[3000:3800,495:510] = True
mask[2008:2080,2614:2691] = True
mask[2500:2535,1045:1070] = True
mask[1073:1081,1892:1904] = True
mask[980:998,1780:1803] = True
mask[3550:3750,3940:3970] = True
mask[2237:2241,3798:3804] = True
mask[1696:1845,2360:2430] = False
mask[1735:1820,2405:2450] = True
mask[1890:1940,2510:2560] = False
mask[1655:1685,2625:2655] = False

#afbeeling wordt getoond en alle sterren worden gezocht. Hieruit wordt dan data extract
plt.imshow(section1, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 
tbl = find_peaks(section1, median+5*std, box_size=11, mask=mask, border_width=10, wcs=wcs)
tbl['peak_value'].info.format = '%.8g'

tbl['id'] = np.arange(len(tbl))
tbl['img'] = '07'
positions = np.transpose((tbl['x_peak'], tbl['y_peak']))
apertures = CircularAperture(positions, r=25)
apertures.plot(color='red', lw=2, alpha=0.5)

#Dit deel voegt de datapunten toe van sterren die uitgemask waren in het eerste gedeelte
plt.subplot(212) 
plt.imshow(section1, origin='lower', cmap='Greys', vmin=median, vmax=median+5*std, interpolation='nearest') 
mask1_2 = np.zeros(section1.shape, dtype=bool)
mask1_2[0:4000,0:4200] = True
mask1_2[1600:2400,2250:4200] = False
tbl1_2 = find_peaks(section1, median+5*std, box_size=50, border_width=10, wcs=wcs, npeaks=4, mask=mask1_2)
positions1_2 = np.transpose((tbl1_2['x_peak'], tbl1_2['y_peak']))
apertures1_2 = CircularAperture(positions1_2, r=25)
apertures1_2.plot(color='red', lw=2, alpha=0.5)

#data formatten en combineren in een enkele tabel
tbl1_2['id'] = np.arange(len(tbl1_2))+len(tbl)
tbl1_2['img'] = '07'
for id in range(len(tbl1_2)):
    tbl.add_row(tbl1_2[id])

overbelichte_sterren = [468, 342, 218, 106, 30, 540, 542, 539, 471, 541]
tbl['OB'] = np.isin(tbl['id'], overbelichte_sterren).astype(int)
print(tbl)

for row in tbl:
    plt.text(int(row['x_peak']), int(row['y_peak']), int(row['id']), color='r', fontsize='medium')

plt.show()



#csv file maken, alleen uncommenten als het nodig is :)
#tbl.write('stars_data_20240307.ecsv', overwrite=True) #file maken met data
