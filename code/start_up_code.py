#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:39:57 2024

@author: heikamp


start up script for astronomy image analysis with python

"""
#import packages

from astropy.io import fits #to read in FITS files
import os # os.path to manipulate file paths 
import glob # finding pathnames (to search for certain fits files in folders and subfolders)
import numpy as np # math applied to arrays (important, no need to read pixel for pixel!)
from matplotlib import pyplot as plt #plot library 
from astropy.visualization import ZScaleInterval #create minimum and maximum Z values for plotting 


#import fits file
data_im ='/Users/ttbjo/programming dings/eerstejaars_project_Praesape/G8_OpenClusters/science_frames/Master_Light_G60s_20240307.fit'

hdu = fits.open(data_im)[0]
header = hdu.header 
data = hdu.data

print(header)

z = ZScaleInterval()
z1,z2 = z.get_limits(data)
plt.imshow(data, vmin=z1, vmax=z2, cmap='gray')


#import list of fits files:
# set path to calibration folder (path_in) and path to save files (path_out)
path_in = " /Users/Documents/DATA/"
path_out=" /Users/Documents/DATA/python/"

file_list = glob.glob(path_in+'/**/'+str("*.fit*"),recursive=True) #/**/ means to search subfolders if recursive = True
sorted_list = sorted(file_list, key=os.path.getmtime) 

hdu_list = fits.open(sorted_list[10])[0] # 10 = 10th item in list
header = hdu_list.header 
data = hdu_list.data


#plate solve image: upload image to nova.astrometry.net/upload


#simple manipulation and save fits file
data = data.T 
data_sub = data_frame1 - data_frame2
data_median = np.median(data)


hdu = fits.PrimaryHDU(data_sub, header=header)
hdu.writeto(path_out + "filename.fits", overwrite=True)



