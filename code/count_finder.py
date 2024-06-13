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
from photutils.aperture import CircularAperture
from photutils.detection import find_peaks

#Cropping area
min_x = 35
max_x = 4500
min_y = 0
max_y = 3900

#Path to dataset containing the pixel coords of the stars
position_data_path = "code\\stars.ecsv"
fit_data_name = "Master_Light_I60s_20240307.fit"


contrast_range = [150.99, 202.11]

#Pixels to add on to the masking radius.
extra_cover = 10

#manually set the radii of certain cases
#id : radius
set_radii = {
    177 : 8,
    178 : 18
}

#Reads a circular area from the image
def read_counts(radius, centroid_x, centroid_y, image, noise):
    # X and Y are flipped for some reason
    area = image[np.clip(int(centroid_y-radius),0,image.shape[0]):np.clip(int(centroid_y+radius+1),0,image.shape[0]),
                 np.clip(int(centroid_x-radius),0,image.shape[1]):np.clip(int(centroid_x+radius+1),0,image.shape[1])]
    pixel_count = 0
    sum = 0
    overexposed = False
    edge_case = False

    if area.shape[1] != radius*2+1 or area.shape[0] != radius*2+1:
        edge_case = True

    for x in range(0,area.shape[1]):
        for y in range(0,area.shape[0]):
            if (x-radius)*(x-radius) + (y-radius)*(y-radius) <= radius*radius:
                sum += area[y,x] - noise
                if area[y,x] >= 65536:
                    overexposed = True
                pixel_count += 1
            else:
                continue
    
    average = sum/pixel_count

    return sum, average, overexposed, edge_case

def cut_background(radius, centroid_x, centroid_y, image):
    # X and Y are flipped for some reason
    image[centroid_y,centroid_x] = -1
    for x in range(np.clip(int(centroid_x-radius),0,image.shape[1]),np.clip(int(centroid_x+radius+1),0,image.shape[1])):
        for y in range(np.clip(int(centroid_y-radius),0,image.shape[0]),np.clip(int(centroid_y+radius+1),0,image.shape[0])):
            #print(x,y)
            #image[y,x] = -1
            if (x-centroid_x)*(x-centroid_x) + (y-centroid_y)*(y-centroid_y) <= radius*radius:
                image[y,x] = -1
            else:
                continue

    return image

def update_mask(radius, centroid_x, centroid_y, mask):
    # X and Y are flipped for some reason
    mask[centroid_y,centroid_x] = -1
    for x in range(np.clip(int(centroid_x-radius),0,mask.shape[1]),np.clip(int(centroid_x+radius+1),0,mask.shape[1])):
        for y in range(np.clip(int(centroid_y-radius),0,mask.shape[0]),np.clip(int(centroid_y+radius+1),0,mask.shape[0])):
            if (x-centroid_x)*(x-centroid_x) + (y-centroid_y)*(y-centroid_y) <= radius*radius:
                mask[y,x] = 1
            else:
                continue

    return mask

def process_background(image, bins):
    histogram, bin_edges = np.histogram(image, bins=bins, range=(0, 65536))
    plt.subplot(111)
    plt.xlim([0,65536])
    plt.ylim(1000)
    plt.plot(bin_edges[0:-1], histogram)

    plt.show()

    median = 0

    for i, j in enumerate(histogram):
        if j == np.max(histogram):
            median = i*(65536/bins)

    mean = sum([i * (65536 / bins) * histogram[i] / (image.size) for i in range(bins)])
    std = sum([(i * (65536 / bins)-mean)**2 * histogram[i] / (image.size) for i in range(bins)])

    return mean, median, std

def get_expected_noise(image):
    value = process_background(image,65536)[1]
    return value

loc = Path(__file__).resolve().parent.parent
loc = PurePath(loc,position_data_path)
print(loc)
star_pos = pd.read_csv(loc, sep=' ', comment='#')

image_file = get_pkg_data_filename(fit_data_name)
hdu = fits.open(image_file)[0]
section = hdu.data[min_y:max_y,min_x:max_x]
background = section
background_mask = np.zeros((section.shape[0],section.shape[1]))

#Noise to be expected for deterimining aperture sizes
expected_noise = process_background(section, 65536)[0]
print(expected_noise)

plt.subplot(111)
plt.imshow(section, origin='lower', cmap='Greys', vmin=contrast_range[0], vmax=contrast_range[1], interpolation='nearest') 

radius = []

for index, row in star_pos.iterrows():
    # X and Y are flipped for some reason
    prev_result = -1000

    for i in range(100):
        result = read_counts(i,int(row['xcentroid']),int(row['ycentroid']),section,expected_noise)
        if prev_result > result[0] or i == 99 or int(row['id']) in set_radii.keys():
            if int(row['id']) in set_radii.keys():
                i = set_radii.get(int(row['id']))
                result = read_counts(i,int(row['xcentroid']),int(row['ycentroid']),section,expected_noise)
            col = 'b'
            if result[2]:
                col = 'r'
            elif result[3]:
                col = 'y'
            plt.gca().add_patch(plt.Circle((row['xcentroid'], row['ycentroid']), i, color=col, lw=2, alpha=0.5, fill=False))
            plt.gca().add_patch(plt.Circle((row['xcentroid'], row['ycentroid']), i*2, color=col, linestyle='dashed', lw=1, alpha=0.5, fill=False))
            plt.gca().add_patch(plt.Circle((row['xcentroid'], row['ycentroid']), i*3, color=col, linestyle='dashed', lw=1, alpha=0.5, fill=False))
            plt.text(row['xcentroid'], row['ycentroid'], int(row['id']), color=col, fontsize='medium')
            radius.append(i)
            #WIP, do actual count using photutils

            background = cut_background(i+extra_cover,int(row['xcentroid']),int(row['ycentroid']),background)
            background_mask = update_mask(i+extra_cover,int(row['xcentroid']),int(row['ycentroid']),background_mask)
            break
        prev_result = result[0]


#plt.subplot(132)
#plt.imshow(background_mask, origin='lower', cmap='Greys', vmin=0, vmax=1, interpolation='nearest') 

#plt.subplot(133)
#plt.imshow(background, origin='lower', cmap='Greys', vmin=contrast_range[0], vmax=contrast_range[1], interpolation='nearest') 

plt.show()

#Get the background noise and standard deviation.
print(process_background(background, 65536))