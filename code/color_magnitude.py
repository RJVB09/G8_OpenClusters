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
import scipy.stats as sps

magnitude_data_output_path_g = "code\\magnitude_data\\magnitudesG.csv"
magnitude_data_output_path_r = "code\\magnitude_data\\magnitudesR.csv"
magnitude_data_output_path_i = "code\\magnitude_data\\magnitudesI.csv"
magnitude_data_output_path_gc = "code\\magnitude_data\\magnitudesClusterNewG.csv"
magnitude_data_output_path_rc = "code\\magnitude_data\\magnitudesClusterNewR.csv"
magnitude_data_output_path_ic = "code\\magnitude_data\\magnitudesClusterNewI.csv"

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

stars_in_both_frames = [351, 317, 493, 516, 112, 153, 9, 27]

test_G = data_GC[(data_GC['star_id'].isin(stars_in_both_frames))]
test_I = data_IC[(data_IC['star_id'].isin(stars_in_both_frames))]
#M = m + 5*(log10(p)+1)

mag_test_G = test_G['magnitudes'] - (5 * (np.log10(1/(test_G['plx']/1000))-1))
mag_test_I = test_I['magnitudes'] - (5 * (np.log10(1/(test_I['plx']/1000))-1))


#I correctiewaarde vinden
G_I_lijst = list(mag_test_G - mag_test_I)
diff_G_I_lijst = []
for n in range(4):
    diff_G_I_lijst.append(G_I_lijst[n] - G_I_lijst[n+4])
correctie_I = np.mean(diff_G_I_lijst)
correctie_I_std = np.std(diff_G_I_lijst)

# als 'img' == 6.0, -corr op mag_IC
data_IC['corr'] = np.where(data_IC['img'] == 7.0, correctie_I, 0)

BS_G = data_GC[data_GC['star_id'] == 539]
BS_I = data_IC[data_IC['star_id'] == 539] 

mag_G_BS = BS_G['magnitudes'] - (5 * (np.log10(1/(BS_G['plx']/1000))-1))
mag_I_BS = BS_I['magnitudes'] - (5 * (np.log10(1/(BS_I['plx']/1000))-1)) + BS_I['corr']

mag_IC = data_IC['magnitudes'] - (5 * (np.log10(1/(data_IC['plx']/1000))-1)) + data_IC['corr']
mag_IC_gcorr = data_IC['magnitudes'] - (5 * (np.log10(1/(data_IC['plx']/1000))-1))
mag_RC = data_RC['magnitudes'] - (5 * (np.log10(1/(data_RC['plx']/1000))-1))
mag_GC = data_GC['magnitudes'] - (5 * (np.log10(1/(data_GC['plx']/1000))-1))

test_I_corr = data_IC[(data_IC['star_id'].isin(stars_in_both_frames))]
mag_test_I_corr = test_I_corr['magnitudes'] - (5 * (np.log10(1/(test_I_corr['plx']/1000))-1)) + test_I_corr['corr']


# isochroon
iso_loc = PurePath(file_loc,"code\\isochroon_0.333ext.dat")
data_iso = pd.read_csv(iso_loc, delim_whitespace=True, comment='#')
data_iso_8 = data_iso[data_iso['logAge'] == 8] 
data_iso_9 = data_iso[data_iso['logAge'] == 9] 
data_iso_corr = data_iso[data_iso['logAge'] == 8.87506] 

iso_mag_G_8 = data_iso_8['gmag']
iso_mag_I_8 = data_iso_8['imag']
iso_mag_G_9 = data_iso_9['gmag']
iso_mag_I_9 = data_iso_9['imag']
iso_mag_G_corr = data_iso_corr['gmag']
iso_mag_I_corr = data_iso_corr['imag']

x = (mag_GC - mag_IC).to_numpy()
x_err = 'ENTER'
y = mag_GC.to_numpy()
y_err = 'ENTER'
isoage = data_iso['logAge'].to_numpy()
isoext = np.full(len(isoage), 0.333)
isoxs = (data_iso['gmag'] - data_iso['imag']).to_numpy()
isoys = data_iso['gmag'].to_numpy()

#SOURCE MITCHEL STOOP
def find_age_with_extinction(x, x_err, y, y_err, isoage, isoext, isoxs, isoys):
    '''
    Find the mean, modal age of individual stars, and the age-prob. distribution with unknown extinction
    This function is a simplified implemenation of the G-function in Jorgensen & Lindegren 2005
    Inspired by D. Guo
    Input:
        x: The x-axis data (np array)
        x_err: The x-axis error on the data (np array)
        y: The y-axis data (np array)
        y_err: The y-axis error on the data (np array)
        isoage: The age column of the isochrone data table (np array)
        isoext:  the extinction column of the isochrone data table (np array)
        isoxs:  The x-axis data for the isochrones (np array)
        isoys:  The y-axis data for the isochrones (np array)

    Return:
        p_all: 2d array containing the best likelihood for the given extinctions and ages for the isochrones
            structured as [extinction i][age j]
            The closer the likelihood to 0, the better. All likelihoods should be negative
    '''
    
    age_list = np.unique(isoage) # all possible ages in a list
    ext_list = np.unique(isoext)

    # make a matrix for the probabilities
    p_all = np.zeros((len(ext_list), len(age_list)))
    
    # iterate over each extinction and age for the isochrones
    for i in range(len(ext_list)):
        av = ext_list[i]
        for j in range(len(age_list)):
            age = age_list[j]
            
            # the x and y normal functions
            # we use vstack on the data so we can make a 2d matrix for each isochrone
            normal_x = sps.norm(np.vstack(x), np.vstack(x_err))
            normal_y = sps.norm(np.vstack(y), np.vstack(y_err))
            
            # calculate the probability for each star for the given isochrone
            # we assume that the covariance matrix is diagonal so we can just multiply the x and y gaussians
            # structured as [star n][isochrone q]
            # amplify by 1e100 so they do not go to infinity when taking the log
            iso_cur_xs = isoxs[(isoage == age)&(isoext == av)]
            iso_cur_ys = isoys[(isoage == age)&(isoext == av)]
            iso_cur_xs = iso_cur_xs[(iso_cur_ys < np.amax(y))&(iso_cur_ys > np.amin(y))]
            iso_cur_ys = iso_cur_ys[(iso_cur_ys < np.amax(y))&(iso_cur_ys > np.amin(y))]
            p_cur = normal_x.pdf(iso_cur_xs)*normal_y.pdf(iso_cur_ys)*1e100

            # sum the probabilities for each star for the given isochrone
            p_cur_sum = np.sum(p_cur, axis=1)
            
            # take the logarithm to avoid large powers
            p_cur_log = np.log10(p_cur_sum)
            
            # set the infinity to -infinity as we want to take the maximum next
            p_cur_log[np.isinf(p_cur_log)] = 0
            
            # sum the maximum probabilities for all stars
            p_sum_log = np.sum(p_cur_log)
            p_all[i,j] = p_sum_log
        
    return p_all



plt.subplot(111)
plt.title('Color-magnitude diagram M44')
#plt.scatter(mag_G - mag_I, mag_G, s = 7, alpha = 1, color = '#7570b3', label='Non-cluster stars')
plt.scatter(iso_mag_G_8-iso_mag_I_8, iso_mag_G_8, s=10, alpha = 0.5, c='green', marker='o', label='isochrone 100 MY')
plt.scatter(iso_mag_G_9-iso_mag_I_9, iso_mag_G_9, s=10, alpha = 0.5, c='blue', marker='o', label='isochrone 1000 MY')
plt.scatter(iso_mag_G_corr-iso_mag_I_corr, iso_mag_G_corr, s=10, alpha = 0.5, c='purple', marker='o', label='isochrone ~750 MY (expection)')
plt.scatter(mag_GC - mag_IC, mag_GC, s = 50, alpha = 1, color = '#d95f02', marker='*',label='Stars of M44')
plt.scatter(mag_G_BS - mag_I_BS, mag_G_BS, s=80, c='blue', marker='*', label='BLUE STRAGGLER JIPPIE')
plt.scatter(mag_test_G - mag_test_I_corr, mag_test_G, s=80, c='black', marker='*', label='DUPLICATE JIPPIE')
plt.legend()
plt.gca().invert_yaxis()
plt.xlabel("G - I")
plt.ylabel("G")
plt.show()


plt.title('Color-color diagram M44')
#plt.scatter(mag_G - mag_I, mag_G, s = 10, alpha = 1, color = '#7570b3', label='Background stars')
plt.scatter(mag_GC - mag_IC_gcorr, mag_GC - mag_RC, s = 10, alpha = 1, color = '#d95f02', label='Stars of M44')
plt.legend()
plt.gca().invert_yaxis()
plt.xlabel("g - i")
plt.ylabel("g - r")
plt.show()

plt.title('Color-color diagram M44')
#plt.scatter(mag_G - mag_I, mag_G, s = 10, alpha = 1, color = '#7570b3', label='Background stars')
plt.scatter(0.98*(mag_GC - mag_RC) + 0.22, 1.09*(mag_RC - mag_IC_gcorr) + 0.22, s = 10, alpha = 1, color = '#d95f02', label='Stars of M44')
plt.legend()
plt.xlabel("B - V")
plt.ylabel("V - R")
plt.show()

