import pandas as pd
import csv

#zet de locatie van hier:
loc = "C:/Users/ttbjo/programming dings/eerstejaars_project_Praesape/G8_OpenClusters/data/asu.tsv"

#data wordt geimporteerd als pandas Dataframe (https://www.geeksforgeeks.org/python-pandas-dataframe/)
data = pd.read_csv(loc, sep=';', comment='#')
print(data)

#als we een bepaalde ster zoeken uit onze data
RA_low = 0.0
RA_high = 360.0
DE_low = 0.0
DE_high = 360.0

#elk datapunt wordt geprint
for index, row in data.iterrows():
    #index, gaia_id, Right Asc, Declination, Parallax, proper motion RA, proper motion 
    print(index, row['GaiaDR3'], row['RA_ICRS'], row['DE_ICRS'], row['Plx'], row['pmRA'], row['pmDE'])
    if index > 1 and float(row['RA_ICRS']) > RA_low and float(row['RA_ICRS']) < RA_high and float(row['DE_ICRS']) > DE_low and float(row['DE_ICRS']) < DE_high:
        print('BINGO') 


