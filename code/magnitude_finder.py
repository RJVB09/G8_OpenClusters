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

#Zero point star, zero_point_counts : counts of the zero point star in the image. zero_point_offset :
#zero_point_counts = 10000
#zero_point_offset = 1

#contrast_range = [150.99, 202.11]
contrast_range = [127.9,217.24]

#Aperture radius
base_radius = 15

manual_zeropoint = True


##calib_parallax = 5.3669 #milli arcseconds
#calib_parallax_err = 0.02
#https://simbad.cds.unistra.fr/simbad/sim-id?Ident=%401108429&Name=BD%2b20%20%202151&submit=submit


def get_magnitude(image_id, overexposed, zeropoint, zeropoint_err, calib_id, calib_mag, calib_mag_err, manual_zeropoint, fit, input_stars):

    file_loc = Path(__file__).resolve().parent.parent
    loc = PurePath(file_loc,input_stars)

    star_pos = pd.read_csv(loc, sep=',', comment='#')
    star_pos = star_pos[star_pos['img'] == image_id]
    star_pos = star_pos[(star_pos['ob'] == overexposed) | (star_pos['id'].isin(calib_id))]

    image_file = PurePath(file_loc, f"science_frames\\aligned\\{fit}")
    hdu = fits.open(image_file)[0]

    section = hdu.data[0:4096,0:4096]

    if image_id == 7.0:
        section = hdu.data[min_y:max_y,min_x:max_x]

    plt.subplot(111)

    positions = np.transpose((star_pos['x_peak'], star_pos['y_peak']))

    for index, row in star_pos.iterrows():
        plt.text(row['x_peak'], row['y_peak'], int(row['id']), color='r', fontsize='medium')


    apertures = CircularAperture(positions, r=base_radius)
    apertures.plot(color='b', lw=2, alpha=0.5)

    annulus_aperture = CircularAnnulus(positions, r_in=base_radius*2,r_out=base_radius*3)
    annulus_aperture.plot(color='r', lw=2, alpha=0.5)

    aperstats = ApertureStats(section, annulus_aperture)
    bkg_mean = aperstats.mean
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

    magnitudes = - 2.5 * np.log10(star_data['corrected_sum'])
    magnitudes_err = 1 / (star_data['corrected_sum'] ** 2) * star_data['corrected_sum_err']

    star_data['inst_magnitudes'] = magnitudes
    star_data['inst_magnitudes_err'] = magnitudes_err
    star_data['star_id'] = star_pos['id']
    if 'plx' in star_pos:
        star_data['plx'] = star_pos['plx']

    star_data_df = star_data.to_pandas()

    #print(calib_mag_data)

    if not manual_zeropoint:
        average_zpt = 0
        for i in range(len(calib_id)):
            calib_mag_data = float(star_data_df[star_data_df['star_id'] == calib_id[i]]['inst_magnitudes'])
            calib_mag_data_err = float(star_data_df[star_data_df['star_id'] == calib_id[i]]['inst_magnitudes_err'])
            average_zpt += calib_mag[i] - calib_mag_data
            zeropoint_err += calib_mag_err[i] ** 2 + calib_mag_data_err ** 2

        zeropoint = average_zpt / len(calib_id)
        zeropoint_err = np.sqrt(zeropoint_err) / len(calib_id)

    star_data['magnitudes'] = zeropoint + magnitudes
    star_data['magnitudes_err'] = np.sqrt(magnitudes_err ** 2 + zeropoint_err ** 2)


    #plt.imshow(section, origin='lower', cmap='Greys', vmin=contrast_range[0], vmax=contrast_range[1], interpolation='nearest') 

    #plt.show()

    star_data.pprint()

    star_data_df = star_data.to_pandas()

    print('zeropoint : ', zeropoint)
    print('zeropoint_err : ', zeropoint_err)
    return zeropoint, zeropoint_err, star_data_df

def dataframes_save(dataframes, path):
    combined = pd.concat(dataframes, axis=0)
    file_loc = Path(__file__).resolve().parent.parent
    output_loc = PurePath(file_loc,path)
    combined.to_csv(output_loc)

#https://simbad.cds.unistra.fr/simbad/sim-id?Ident=%401108429&Name=BD%2b20%20%202151&submit=submit
calib_id_7 = [203]
calib_magGRI_7 = [[10.619],[10.032],[9.843]]
calib_magGRI_7_err = [[0.073],[0.054],[0.069]]


#https://simbad.cds.unistra.fr/simbad/sim-id?Ident=%401109963&Name=BD%2b20%20%202170&submit=submit
#https://simbad.cds.unistra.fr/simbad/sim-id?Ident=%4011682158&Name=BD%2b20%20%202160&submit=submit
#https://simbad.cds.unistra.fr/simbad/sim-id?Ident=%401110935&Name=Cl*%20NGC%202632%20%20%20%20%20%20S%20%20%20%20%20106&submit=submit
calib_id_6 = [224,198]
calib_magGRI_6 = [[13.962,11.991],[13.504,11.718],[13.290,11.664]]
calib_magGRI_6_err = [[0.009,0.001],[0.01,0.001],[0.012,0.001]]
"""
calib_id_6 = [281,224,198]
calib_magGRI_6 = [[13.643,13.962,11.991],[13.393,13.504,11.718],[13.184,13.290,11.664]]
calib_magGRI_6_err = [[0.01,0.009,0.001],[0.011,0.01,0.001],[0.012,0.012,0.001]]
"""

mag_g_o0_i7_c0 = get_magnitude(7.0, 0, 0, 0, calib_id_7, calib_magGRI_7[0], calib_magGRI_7_err[0], False, "Master_Light_G60s_20240307.fit", "code\\stars_data_20240307.csv")
mag_r_o0_i7_c0 = get_magnitude(7.0, 0, 0, 0, calib_id_7, calib_magGRI_7[1], calib_magGRI_7_err[1], False, "Master_Light_R60s_20240307.fit", "code\\stars_data_20240307.csv")
mag_i_o0_i7_c0 = get_magnitude(7.0, 0, 0, 0, calib_id_7, calib_magGRI_7[2], calib_magGRI_7_err[2], False, "Master_Light_I60s_20240307.fit", "code\\stars_data_20240307.csv")

mag_g_o1_i7_c0 = get_magnitude(7.0, 1, 0, 0, calib_id_7, calib_magGRI_7[0], calib_magGRI_7_err[0], False, "Master_Light_G4s_20240307.fit", "code\\stars_data_20240307.csv")
mag_r_o1_i7_c0 = get_magnitude(7.0, 1, 0, 0, calib_id_7, calib_magGRI_7[1], calib_magGRI_7_err[1], False, "Master_Light_R4s_20240307.fit", "code\\stars_data_20240307.csv")
mag_i_o1_i7_c0 = get_magnitude(7.0, 1, 0, 0, calib_id_7, calib_magGRI_7[2], calib_magGRI_7_err[2], False, "Master_Light_I5s_20240307.fit", "code\\stars_data_20240307.csv")

mag_g_o0_i6_c0 = get_magnitude(6.0, 0, 0, 0, calib_id_6, calib_magGRI_6[0], calib_magGRI_6_err[0], False, "Master_Light_G60s_20240306.fit", "code\\stars_data_20240306.csv")
mag_r_o0_i6_c0 = get_magnitude(6.0, 0, 0, 0, calib_id_6, calib_magGRI_6[1], calib_magGRI_6_err[1], False, "Master_Light_R60s_20240306.fit", "code\\stars_data_20240306.csv")
mag_i_o0_i6_c0 = get_magnitude(6.0, 0, 0, 0, calib_id_6, calib_magGRI_6[2], calib_magGRI_6_err[2], False, "Master_Light_I60s_20240306.fit", "code\\stars_data_20240306.csv")

mag_g_o1_i6_c0 = get_magnitude(6.0, 1, 0, 0, calib_id_6, calib_magGRI_6[0], calib_magGRI_6_err[0], False, "Master_Light_G5s_20240306.fit", "code\\stars_data_20240306.csv")
mag_r_o1_i6_c0 = get_magnitude(6.0, 1, 0, 0, calib_id_6, calib_magGRI_6[1], calib_magGRI_6_err[1], False, "Master_Light_R4s_20240306.fit", "code\\stars_data_20240306.csv")
mag_i_o1_i6_c0 = get_magnitude(6.0, 1, 0, 0, calib_id_6, calib_magGRI_6[2], calib_magGRI_6_err[2], False, "Master_Light_I5s_20240306.fit", "code\\stars_data_20240306.csv")

dataframes_save([mag_g_o0_i7_c0[2], mag_g_o1_i7_c0[2], mag_g_o0_i6_c0[2], mag_g_o1_i6_c0[2]],"code\\magnitude_data\\magnitudesG.csv")
dataframes_save([mag_r_o0_i7_c0[2], mag_r_o1_i7_c0[2], mag_r_o0_i6_c0[2], mag_r_o1_i6_c0[2]],"code\\magnitude_data\\magnitudesR.csv")
dataframes_save([mag_i_o0_i7_c0[2], mag_i_o1_i7_c0[2], mag_i_o0_i6_c0[2], mag_i_o1_i6_c0[2]],"code\\magnitude_data\\magnitudesI.csv")

#get_magnitude(7.0, 0, 25.212015560856376, 0.0730000000000001, 203, 10.619, 0.073, True, "Master_Light_G60s_20240307.fit", "code\\member_stars.csv", "code\\magnitudeGC.csv")


# CLUSTER STARS

member_stars = "code\\member_stars_test.csv"

mag_g_o0_i7_c1 = get_magnitude(7.0, 0, mag_g_o0_i7_c0[0], mag_g_o0_i7_c0[1], [], [], [], True, "Master_Light_G60s_20240307.fit", member_stars)
mag_r_o0_i7_c1 = get_magnitude(7.0, 0, mag_r_o0_i7_c0[0], mag_r_o0_i7_c0[1], [], [], [], True, "Master_Light_R60s_20240307.fit", member_stars)
mag_i_o0_i7_c1 = get_magnitude(7.0, 0, mag_i_o0_i7_c0[0], mag_i_o0_i7_c0[1], [], [], [], True, "Master_Light_I60s_20240307.fit", member_stars)

mag_g_o1_i7_c1 = get_magnitude(7.0, 1, mag_g_o1_i7_c0[0], mag_g_o1_i7_c0[1], [], [], [], True, "Master_Light_G4s_20240307.fit", member_stars)
mag_r_o1_i7_c1 = get_magnitude(7.0, 1, mag_r_o1_i7_c0[0], mag_r_o1_i7_c0[1], [], [], [], True, "Master_Light_R4s_20240307.fit", member_stars)
mag_i_o1_i7_c1 = get_magnitude(7.0, 1, mag_i_o1_i7_c0[0], mag_i_o1_i7_c0[1], [], [], [], True, "Master_Light_I5s_20240307.fit", member_stars)

mag_g_o0_i6_c1 = get_magnitude(6.0, 0, mag_g_o0_i6_c0[0], mag_g_o0_i6_c0[1], [], [], [], True, "Master_Light_G60s_20240306.fit", member_stars)
mag_r_o0_i6_c1 = get_magnitude(6.0, 0, mag_r_o0_i6_c0[0], mag_r_o0_i6_c0[1], [], [], [], True, "Master_Light_R60s_20240306.fit", member_stars)
mag_i_o0_i6_c1 = get_magnitude(6.0, 0, mag_i_o0_i6_c0[0], mag_i_o0_i6_c0[1], [], [], [], True, "Master_Light_I60s_20240306.fit", member_stars)

mag_g_o1_i6_c1 = get_magnitude(6.0, 1, mag_g_o1_i6_c0[0], mag_g_o1_i6_c0[1], [], [], [], True, "Master_Light_G5s_20240306.fit", member_stars)
mag_r_o1_i6_c1 = get_magnitude(6.0, 1, mag_r_o1_i6_c0[0], mag_r_o1_i6_c0[1], [], [], [], True, "Master_Light_R4s_20240306.fit", member_stars)
mag_i_o1_i6_c1 = get_magnitude(6.0, 1, mag_i_o1_i6_c0[0], mag_i_o1_i6_c0[1], [], [], [], True, "Master_Light_I5s_20240306.fit", member_stars)

dataframes_save([mag_g_o0_i7_c1[2], mag_g_o1_i7_c1[2], mag_g_o0_i6_c1[2], mag_g_o1_i6_c1[2]],"code\\magnitude_data\\magnitudesClusterNewG.csv")
dataframes_save([mag_r_o0_i7_c1[2], mag_r_o1_i7_c1[2], mag_r_o0_i6_c1[2], mag_r_o1_i6_c1[2]],"code\\magnitude_data\\magnitudesClusterNewR.csv")
dataframes_save([mag_i_o0_i7_c1[2], mag_i_o1_i7_c1[2], mag_i_o0_i6_c1[2], mag_i_o1_i6_c1[2]],"code\\magnitude_data\\magnitudesClusterNewI.csv")