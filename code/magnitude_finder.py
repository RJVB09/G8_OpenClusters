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

loc = Path(__file__).resolve().parent.parent
loc = PurePath(loc,"data\\asu.tsv")
print(loc)

star_pos = pd.read_csv(loc, sep=';', comment='#')
print(star_pos)

for index, row in star_pos.iterrows():
    print(row['c1'], row['c2'])