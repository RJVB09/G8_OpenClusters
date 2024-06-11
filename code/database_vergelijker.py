import pandas as pd
from pathlib import Path
from pathlib import PurePath

#flexible file location
main_loc = Path(__file__).resolve().parent.parent
loc = PurePath(main_loc,"data\\asu.tsv")
print(loc)

#data wordt geimporteerd als pandas Dataframe (https://www.geeksforgeeks.org/python-pandas-dataframe/)
data_db = pd.read_csv(loc, sep=';', comment='#')
print(data_db)

#als we een bepaalde ster zoeken uit onze data
RA_low = 0.0
RA_high = 360.0
DE_low = 0.0
DE_high = 360.0

#elk datapunt wordt geprint
for index, row in data_db.iterrows():
    #index, gaia_id, Right Asc, Declination, Parallax, proper motion RA, proper motion 
    print(index, row['GaiaDR3'], row['RA_ICRS'], row['DE_ICRS'], row['Plx'], row['pmRA'], row['pmDE'])
    #if index > 1 and float(row['RA_ICRS']) > RA_low and float(row['RA_ICRS']) < RA_high and float(row['DE_ICRS']) > DE_low and float(row['DE_ICRS']) < DE_high:
        #print('BINGO') 

#import hier database met meer sterren van Mitchel, dan kunnen we data dubbelchecken en vergelijken

#import hier de twee sterdatabases die we zelf hebben!
loc_stars_06 = PurePath(main_loc,"code\\stars_data_20240306.ecsv")
loc_stars_07 = PurePath(main_loc,"code\\stars_data_20240307.ecsv")
data_stars_06 = pd.read_csv(loc_stars_06, sep=';', comment='#')
data_stars_07 = pd.read_csv(loc_stars_07, sep=';', comment='#')

print(data_stars_06)
for index, row in data_stars_06.iterrows():
    print(row)
    #print(row['id'], row['skycoord_peak.ra'], row['skycoord_peak.dec'])