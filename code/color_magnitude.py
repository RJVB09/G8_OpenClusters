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

magnitude_data_output_path_g = "code\\magnitudeG.csv"
magnitude_data_output_path_r = "code\\magnitudeR.csv"
magnitude_data_output_path_i = "code\\magnitudeI.csv"
magnitude_data_output_path_gc = "code\\magnitudeGC.csv"
magnitude_data_output_path_rc = "code\\magnitudeRC.csv"
magnitude_data_output_path_ic = "code\\magnitudeIC.csv"

file_loc = Path(__file__).resolve().parent.parent
locG = PurePath(file_loc,magnitude_data_output_path_g)
locR = PurePath(file_loc,magnitude_data_output_path_r)
locI = PurePath(file_loc,magnitude_data_output_path_i)
locGC = PurePath(file_loc,magnitude_data_output_path_gc)
locRC = PurePath(file_loc,magnitude_data_output_path_rc)
locIC = PurePath(file_loc,magnitude_data_output_path_ic)

data_G = pd.read_csv(locG, sep=',', comment='#')
data_R = pd.read_csv(locR, sep=',', comment='#')
data_I = pd.read_csv(locI, sep=',', comment='#')
data_GC = pd.read_csv(locGC, sep=',', comment='#')
data_RC = pd.read_csv(locRC, sep=',', comment='#')
data_IC = pd.read_csv(locIC, sep=',', comment='#')

mag_G = data_G['magnitudes']
mag_R = data_R['magnitudes']
mag_I = data_I['magnitudes']
mag_IC = data_IC['magnitudes']
mag_RC = data_RC['magnitudes']
mag_GC = data_GC['magnitudes']


# isochroon
iso_loc = PurePath(file_loc,"code\\isochroon_test.dat")
data_iso = pd.read_csv(iso_loc, delim_whitespace=True, comment='#')
print(data_iso)
print(data_iso.columns)
data_iso_8 = data_iso[data_iso['logAge'] == 8] 
data_iso_9 = data_iso[data_iso['logAge'] == 9] 
data_iso_corr = data_iso[data_iso['logAge'] == 8.87506] 
mag_corr = 6
iso_mag_G_8 = data_iso_8['gmag'] +mag_corr
iso_mag_I_8 = data_iso_8['imag'] +mag_corr
iso_mag_G_9 = data_iso_9['gmag'] +mag_corr
iso_mag_I_9 = data_iso_9['imag'] +mag_corr
iso_mag_G_corr = data_iso_corr['gmag'] +mag_corr
iso_mag_I_corr = data_iso_corr['imag'] +mag_corr



plt.subplot(111)
plt.title('Color-magnitude diagram M44')
#plt.scatter(mag_G - mag_I, mag_G, s = 7, alpha = 1, color = '#7570b3', label='Non-cluster stars')
plt.scatter(iso_mag_G_8-iso_mag_I_8, iso_mag_I_8, s=3, alpha = 0.3, c='green', marker='*')
plt.scatter(iso_mag_G_9-iso_mag_I_9, iso_mag_I_9, s=3, alpha = 0.3, c='blue', marker='*')
plt.scatter(iso_mag_G_corr-iso_mag_I_corr, iso_mag_I_corr, s=3, alpha = 0.3, c='purple', marker='*')
plt.scatter(mag_GC - mag_IC, mag_GC, s = 10, alpha = 1, color = '#d95f02', label='Stars of M44')
plt.legend()
plt.gca().invert_yaxis()
plt.xlabel("G - I")
plt.ylabel("G")
plt.show()

