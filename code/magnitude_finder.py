import pandas as pd
from pathlib import Path
from pathlib import PurePath
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
from matplotlib.colors import LogNorm
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats,aperture_photometry
from photutils.detection import find_peaks
from acstools import acszpt

#Cropping area
min_x = 35
max_x = 4500
min_y = 0
max_y = 3900

#Path to dataset containing the pixel coords of the stars
position_data_path = "code\\small_star_dataset.csv"
magnitude_data_output_path = "code\\magnitudes.csv"
fit_data_name = "Master_Light_I60s_20240307.fit"

#Zero point star, zero_point_counts : counts of the zero point star in the image. zero_point_offset :
#zero_point_counts = 10000
#zero_point_offset = 1

#Header data
date = '2023-03-07'
instrument = 'WFC' #!
filter = 'F775W' #!
exp_time = 60

contrast_range = [150.99, 202.11]
#contrast_range = [99.073,127.31]

#Aperture radius
base_radius = 15

#calibration

file_loc = Path(__file__).resolve().parent.parent
loc = PurePath(file_loc,position_data_path)
print(loc)
output_loc = PurePath(file_loc,magnitude_data_output_path)

star_pos = pd.read_csv(loc, sep=' ', comment='#')

image_file = get_pkg_data_filename(fit_data_name)
hdu = fits.open(image_file)[0]
section = hdu.data[min_y:max_y,min_x:max_x]

plt.subplot(111)

positions = np.transpose((star_pos['xcentroid'], star_pos['ycentroid']))

for index, row in star_pos.iterrows():
    # X and Y are flipped for some reason
    prev_result = -1000

    plt.text(row['xcentroid'], row['ycentroid'], int(row['id']), color='r', fontsize='medium')


apertures = CircularAperture(positions, r=base_radius)
apertures.plot(color='b', lw=2, alpha=0.5)

annulus_aperture = CircularAnnulus(positions, r_in=base_radius*2,r_out=base_radius*3)
annulus_aperture.plot(color='r', lw=2, alpha=0.5)

aperstats = ApertureStats(section, annulus_aperture)
bkg_mean = aperstats.mean
#bkg_err = np.sqrt(aperstats.mean) #!
bkg_err = aperstats.std #!
aperture_area = apertures.area_overlap(section)
total_bkg = bkg_mean*aperture_area
total_bkg_err = bkg_err*aperture_area

star_data = aperture_photometry(section, apertures, np.sqrt(section))

star_data['total_bkg'] = total_bkg
star_data['total_bkg_err'] = total_bkg_err
star_data['corrected_sum'] = star_data['aperture_sum'] - total_bkg
star_data['corrected_sum_err'] = np.sqrt(total_bkg_err ** 2 + star_data['aperture_sum_err'] ** 2)

for col in star_data.colnames:
    star_data[col].info.format = '%.8g'

#print(star_data)

q = acszpt.Query(date = date, detector = instrument)
zpt_table = q.fetch()
q_filter = acszpt.Query(date = date, detector = instrument, filt = filter)
filter_zpt = q_filter.fetch()

print(filter_zpt)
print(np.float16(filter_zpt['VEGAmag'][0]))

zeropoint = np.float16(filter_zpt['VEGAmag'][0])

magnitudes = zeropoint - 2.5 * np.log10(star_data['corrected_sum'] / exp_time)
magnitudes_err = 1 / (star_data['corrected_sum'] ** 2) * star_data['corrected_sum_err']

star_data['abs_magnitudes'] = magnitudes
star_data['abs_magnitudes_err'] = magnitudes_err

star_data.pprint()
star_data.to_pandas().to_csv(output_loc)

plt.imshow(section, origin='lower', cmap='Greys', vmin=contrast_range[0], vmax=contrast_range[1], interpolation='nearest') 

plt.show()


