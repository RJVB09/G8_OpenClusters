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

#Reads a circular area from the image
def read_counts(radius, centroid_x, centroid_y, image, noise):
    # X and Y are flipped for some reason
    area = image[int(centroid_y-radius):int(centroid_y+radius+1),int(centroid_x-radius):int(centroid_x+radius+1)]
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



loc = Path(__file__).resolve().parent.parent
loc = PurePath(loc,"code\\stars.ecsv")
print(loc)

star_pos = pd.read_csv(loc, sep=' ', comment='#')
#print(star_pos)

image_file = get_pkg_data_filename('20240307_data.fits')
hdu = fits.open(image_file)[0]
section1 = hdu.data[0:3900,35:4500]
#positions = np.transpose((star_pos['xcentroid'], star_pos['ycentroid']))
#apertures = CircularAperture(positions, r=25)
#apertures.plot(color='red', lw=2, alpha=0.5)

#print(section1)

plt.subplot(111)
plt.imshow(section1, origin='lower', cmap='Greys', vmin=512, vmax=674, interpolation='nearest') 

positions = np.transpose((star_pos['xcentroid'], star_pos['ycentroid']))


for index, row in star_pos.iterrows():
    # X and Y are flipped for some reason
    prev_result = -1000
    for i in range(100):
        result = read_counts(i,int(row['xcentroid']),int(row['ycentroid']),section1,544)[0]
        if prev_result > result or i == 99:
            
            plt.gca().add_patch(plt.Circle((row['xcentroid'],row['ycentroid']),i,color='b', lw=2, alpha=0.5, fill=False))
            print(result, i, row['id'])
            break
        prev_result = result

plt.show()

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


print(process_background(section1,65536))