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

magnitude_data_output_path_g = "code\\magnitude_data\\magnitudesG.csv"
magnitude_data_output_path_r = "code\\magnitude_data\\magnitudesR.csv"
magnitude_data_output_path_i = "code\\magnitude_data\\magnitudesI.csv"
magnitude_data_output_path_gc = "code\\magnitude_data\\magnitudesClusterG.csv"
magnitude_data_output_path_rc = "code\\magnitude_data\\magnitudesClusterR.csv"
magnitude_data_output_path_ic = "code\\magnitude_data\\magnitudesClusterI.csv"

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

BS_G = data_GC[data_GC['star_id'] == 539]
BS_I = data_IC[data_IC['star_id'] == 539] 

#M = m + 5*(log10(p)+1)
mag_G_BS = BS_G['magnitudes'] - (5 * (np.log10(1/(BS_G['plx']/1000))-1))
mag_I_BS = BS_I['magnitudes'] - (5 * (np.log10(1/(BS_I['plx']/1000))-1))

#mag_G = data_G['magnitudes'] + 5 * (np.log10(data_G['plx'])+1)
#mag_R = data_R['magnitudes'] + 5 * (np.log10(data_R['plx'])+1)
#mag_I = data_I['magnitudes'] + 5 * (np.log10(data_I['plx'])+1)
mag_IC = data_IC['magnitudes'] - (5 * (np.log10(1/(data_IC['plx']/1000))-1))
mag_RC = data_RC['magnitudes'] - (5 * (np.log10(1/(data_RC['plx']/1000))-1))
mag_GC = data_GC['magnitudes'] - (5 * (np.log10(1/(data_GC['plx']/1000))-1))

# isochroon
iso_loc = PurePath(file_loc,"code\\isochroon_test.dat")
data_iso = pd.read_csv(iso_loc, delim_whitespace=True, comment='#')
print(data_iso)
print(data_iso.columns)
data_iso_8 = data_iso[data_iso['logAge'] == 8] 
data_iso_9 = data_iso[data_iso['logAge'] == 9] 
data_iso_corr = data_iso[data_iso['logAge'] == 8.87506] 

iso_mag_G_8 = data_iso_8['gmag']
iso_mag_I_8 = data_iso_8['imag']
iso_mag_G_9 = data_iso_9['gmag']
iso_mag_I_9 = data_iso_9['imag']
iso_mag_G_corr = data_iso_corr['gmag']
iso_mag_I_corr = data_iso_corr['imag']



plt.subplot(111)
plt.title('Color-magnitude diagram M44')
#plt.scatter(mag_G - mag_I, mag_G, s = 7, alpha = 1, color = '#7570b3', label='Non-cluster stars')
plt.scatter(iso_mag_G_8-iso_mag_I_8, iso_mag_G_8, s=10, alpha = 0.5, c='green', marker='o', label='isochrone 100 MY')
plt.scatter(iso_mag_G_9-iso_mag_I_9, iso_mag_G_9, s=10, alpha = 0.5, c='blue', marker='o', label='isochrone 1000 MY')
plt.scatter(iso_mag_G_corr-iso_mag_I_corr, iso_mag_G_corr, s=10, alpha = 0.5, c='purple', marker='o', label='isochrone ~750 MY (expection)')
plt.scatter(mag_GC - mag_IC, mag_GC, s = 50, alpha = 1, color = '#d95f02', marker='*',label='Stars of M44')
plt.scatter(mag_G_BS - mag_I_BS, mag_G_BS, s=80, c='blue', marker='*', label='BLUE STRAGGLER JIPPIE')
plt.legend()
plt.gca().invert_yaxis()
plt.xlabel("G - I")
plt.ylabel("G")
plt.show()

