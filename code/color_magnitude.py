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

# lijst van de sterren die in beide frames voorkomen
test_G = data_GC[(data_GC['star_id'].isin(stars_in_both_frames))]
test_I = data_IC[(data_IC['star_id'].isin(stars_in_both_frames))]
#M = m + 5*(log10(p)+1)

#absolute magnitudes
parallax_error = 0.03
mag_test_G = test_G['magnitudes'] - (5 * (np.log10(100/(test_G['plx']))))
mag_test_I = test_I['magnitudes'] - (5 * (np.log10(100/(test_I['plx']))))
mag_test_G_err = np.sqrt(test_G['magnitudes_err'] ** 2 + (5 / (np.log(10) * test_G['plx'])) ** 2 * parallax_error ** 2)
mag_test_I_err = np.sqrt(test_G['magnitudes_err'] ** 2 + (5 / (np.log(10) * test_I['plx'])) ** 2 * parallax_error ** 2)
print(test_I)
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
data_IC['corr_err'] = np.where(data_IC['img'] == 7.0, correctie_I_std, 0)
print(data_IC)

BS_G = data_GC[data_GC['star_id'] == 539]
BS_I = data_IC[data_IC['star_id'] == 539] 

mag_G_BS = BS_G['magnitudes'] - (5 * (np.log10(100/BS_G['plx'])))
mag_I_BS = BS_I['magnitudes'] - (5 * (np.log10(100/BS_I['plx']))) + BS_I['corr']

mag_G_BS_err = np.sqrt(BS_G['magnitudes_err'] ** 2 + (5 / (np.log(10) * BS_G['plx'])) ** 2 * parallax_error ** 2)
mag_I_BS_err = np.sqrt(BS_I['magnitudes_err'] ** 2 + (5 / (np.log(10) * BS_I['plx'])) ** 2 * parallax_error ** 2 + BS_I['corr_err'] ** 2)

mag_IC = data_IC['magnitudes'] - (5 * (np.log10(100/data_IC['plx']))) + data_IC['corr']
mag_IC_gcorr = data_IC['magnitudes'] - (5 * (np.log10(100/data_IC['plx'])))

mag_IC_err = np.sqrt(data_IC['magnitudes_err'] ** 2 + (5 / (np.log(10) * data_IC['plx'])) ** 2 * parallax_error ** 2 + data_IC['corr_err'] ** 2)
mag_IC_gcorr_err = np.sqrt(data_IC['magnitudes_err'] ** 2 + (5 / (np.log(10) * data_IC['plx'])) ** 2 * parallax_error ** 2)

mag_RC = data_RC['magnitudes'] - (5 * (np.log10(100/data_RC['plx'])))
mag_GC = data_GC['magnitudes'] - (5 * (np.log10(100/data_GC['plx'])))

mag_RC_err = np.sqrt(data_RC['magnitudes_err'] ** 2 + (5 / (np.log(10) * data_RC['plx'])) ** 2 * parallax_error ** 2)
mag_GC_err = np.sqrt(data_GC['magnitudes_err'] ** 2 + (5 / (np.log(10) * data_GC['plx'])) ** 2 * parallax_error ** 2)

test_I_corr = data_IC[(data_IC['star_id'].isin(stars_in_both_frames))]
mag_test_I_corr = test_I_corr['magnitudes'] - (5 * (np.log10(100/test_I_corr['plx']))) + test_I_corr['corr']

mag_test_I_corr_err = np.sqrt(test_I_corr['magnitudes_err'] ** 2 + (5 / (np.log(10) * test_I_corr['plx'])) ** 2 * parallax_error ** 2 + test_I_corr['corr_err'] ** 2)


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

cluster_stars = pd.DataFrame()
cluster_stars['G-I'] = mag_GC - mag_IC
cluster_stars['G'] = mag_GC
cluster_stars['G-I_err'] = np.sqrt(mag_GC_err ** 2 + mag_IC_err ** 2)
cluster_stars['G_err'] = mag_GC_err

cluster_stars_selected = cluster_stars[cluster_stars['G-I_err'] <= 0.5]
cluster_stars_rest = cluster_stars[cluster_stars['G-I_err'] > 0.5]
x = (mag_GC - mag_IC).to_numpy()
x_err = (np.sqrt(mag_GC_err**2 + mag_IC_err**2)).to_numpy()
y = mag_GC.to_numpy()
y_err = mag_GC_err.to_numpy()
isoage = data_iso['logAge'].to_numpy()
isoext = np.full(len(isoage), 0.333)
isoxs = (data_iso['gmag'] - data_iso['imag']).to_numpy()
isoys = data_iso['gmag'].to_numpy()
nan_ids = np.isnan(x)
print(nan_ids)
x_filt = x[~nan_ids]
x_err_filt = x_err[~nan_ids]
y_filt = y[~nan_ids]
y_err_filt = y_err[~nan_ids]
#too high error
err_ids = np.less_equal(x_filt, np.full(len(x_filt), 0.5))
x_filt_2 = x_filt[~err_ids]
x_err_filt_2 = x_err_filt[~err_ids]
y_filt_2 = y_filt[~err_ids]
y_err_filt_2 = y_err_filt[~err_ids]

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
p_all = find_age_with_extinction(x_filt, x_err_filt, y_filt, y_err_filt, isoage, isoext, isoxs, isoys)
p_all_err_filtered = find_age_with_extinction(x_filt_2, x_err_filt_2, y_filt_2, y_err_filt_2, isoage, isoext, isoxs, isoys)
p_all = p_all - np.amin(p_all)
p_all = p_all / np.amax(p_all)
p_all_err_filtered = p_all_err_filtered - np.amin(p_all_err_filtered)
p_all_err_filtered = p_all_err_filtered / np.amax(p_all_err_filtered)
plt.plot(np.unique(isoage),p_all.reshape(-1), label='Propability all members')
plt.plot(np.unique(isoage),p_all_err_filtered.reshape(-1), label='Propability magnitude error $\leq$0.5')
plt.vlines([np.log10(600000000),np.log10(700000000)],0,1,'black', linestyles='dashed', label='Expectation 600-700 MY')
plt.ylabel('Propability')
plt.xlabel('$log_{10}$(Age) (year)')
plt.legend()
plt.show()

plt.subplot(111)
plt.title('Color-magnitude diagram M44')
#plt.scatter(mag_G - mag_I, mag_G, s = 7, alpha = 1, color = '#7570b3', label='Non-cluster stars')
#plt.scatter(iso_mag_G_8-iso_mag_I_8, iso_mag_G_8, s=10, alpha = 0.5, c='green', marker='o', label='isochrone 100 MY')
#plt.scatter(iso_mag_G_9-iso_mag_I_9, iso_mag_G_9, s=10, alpha = 0.5, c='blue', marker='o', label='isochrone 1000 MY')
plt.scatter(iso_mag_G_corr-iso_mag_I_corr, iso_mag_G_corr, s=10, alpha = 0.5, c='purple', marker='o', label='isochrone ~750 MY')
plt.errorbar(cluster_stars_selected['G-I'], cluster_stars_selected['G'], xerr= cluster_stars_selected['G-I_err'], yerr=cluster_stars_selected['G_err'], color = '#d95f02', fmt="*",label='Stars of M44')
plt.errorbar(cluster_stars_rest['G-I'], cluster_stars_rest['G'], xerr= cluster_stars_rest['G-I_err'], yerr=cluster_stars_rest['G_err'], color = '#666666', fmt="*",label='Stars of M44 (g-i error > 0.5)', alpha = 0.1)
#plt.scatter(mag_GC - mag_IC, mag_GC, s = 50, alpha = 1, color = '#d95f02', marker='*',label='Stars of M44')
plt.errorbar(mag_G_BS - mag_I_BS, mag_G_BS, xerr= np.sqrt(mag_G_BS_err ** 2 + mag_I_BS_err ** 2), yerr=mag_G_BS_err, color = '#7570b3', fmt="*")

plt.errorbar(mag_test_G - mag_test_I_corr, mag_test_G, xerr= np.sqrt(mag_test_G_err ** 2 + mag_test_I_corr_err ** 2), yerr=mag_test_G_err, color = '#666666', fmt="*",label='Stars in both frames')
plt.scatter(mag_G_BS - mag_I_BS, mag_G_BS, s=100, c='#7570b3', marker='*', label='Literature blue straggler')
#plt.scatter(mag_test_G - mag_test_I_corr, mag_test_G, s=80, c='black', marker='*', label='Stars in both frames')
plt.legend()
plt.gca().invert_yaxis()
plt.xlabel("G - I")
plt.ylabel("G")
plt.show()

data_GC['mag_GC - mag_IC'] = mag_GC - mag_IC
data_GC['np.sqrt(mag_GC_err ** 2 + mag_IC_err ** 2)'] = np.sqrt(mag_GC_err ** 2 + mag_IC_err ** 2)
data_GC['mag_GC'] = mag_GC
data_GC['mag_GC_err'] = mag_GC_err


#file_loc = Path(__file__).resolve().parent.parent
#output_loc = PurePath(file_loc,"code\\AAAAA.csv")
#data_GC.to_csv(output_loc)

#print(data_GC)


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

