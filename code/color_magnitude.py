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

magnitude_data_output_path_g = "code\\magnitudesG.csv"
magnitude_data_output_path_r = "code\\magnitudesR.csv"
magnitude_data_output_path_i = "code\\magnitudesI.csv"

file_loc = Path(__file__).resolve().parent.parent
locG = PurePath(file_loc,magnitude_data_output_path_g)
locR = PurePath(file_loc,magnitude_data_output_path_r)
locI = PurePath(file_loc,magnitude_data_output_path_i)

data_G = pd.read_csv(locG, sep=',', comment='#')
data_R = pd.read_csv(locR, sep=',', comment='#')
data_I = pd.read_csv(locI, sep=',', comment='#')

mag_G = data_G['magnitudes']
mag_R = data_R['magnitudes']
mag_I = data_I['magnitudes']

plt.subplot(111)
plt.scatter(mag_G - mag_I, mag_G, s = 5, alpha = 0.5, color = '#d95f02')
plt.gca().invert_yaxis()
plt.xlabel("G - I")
plt.ylabel("G")
plt.show()

