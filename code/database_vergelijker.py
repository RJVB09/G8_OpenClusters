import pandas as pd
from pathlib import Path
from pathlib import PurePath
from IPython.display import display
import matplotlib.pyplot as plt

#flexible file location
main_loc = Path(__file__).resolve().parent.parent
loc = PurePath(main_loc,"data\\asu.tsv")

#data wordt geimporteerd als pandas Dataframe (https://www.geeksforgeeks.org/python-pandas-dataframe/)
data_db = pd.read_csv(loc, sep=';', comment='#')
data_db = data_db.drop([0,1])
data_db['DE_ICRS'] = data_db['DE_ICRS'].str.replace('+', '')
data_db['DE_ICRS'] = data_db['DE_ICRS'].astype(float)
data_db['RA_ICRS'] = data_db['RA_ICRS'].astype(float)

#import hier database met meer sterren van Mitchel, dan kunnen we data dubbelchecken en vergelijken
loc_gaudain = PurePath(main_loc,"data\\cantat_gaudain2018.tsv")
data_goudain = pd.read_csv(loc_gaudain, sep='\t', comment='#')
data_goudain = data_goudain.drop([0,1])
data_goudain['DE_ICRS'] = data_goudain['DE_ICRS'].str.replace('+', '')
data_goudain['DE_ICRS'] = data_goudain['DE_ICRS'].astype(float)
data_goudain['RA_ICRS'] = data_goudain['RA_ICRS'].astype(float)
display(data_goudain)

#import hier de twee sterdatabases die we zelf hebben!
loc_stars_06 = PurePath(main_loc,"code\\stars_data_20240306.ecsv")
loc_stars_07 = PurePath(main_loc,"code\\stars_data_20240307.ecsv")
data_stars_06 = pd.read_csv(loc_stars_06, sep=' ', comment='#')
data_stars_07 = pd.read_csv(loc_stars_07, sep=' ', comment='#')

display(data_stars_06)
display(data_stars_07)
print(data_stars_07['skycoord_peak.ra'])
print(data_db['DE_ICRS'])
print(data_db['RA_ICRS'])

plt.scatter(data_stars_06['skycoord_peak.ra'],data_stars_06['skycoord_peak.dec'],color='red')
plt.scatter(data_stars_07['skycoord_peak.ra'],data_stars_07['skycoord_peak.dec'],color='blue')
plt.scatter(x=data_db['RA_ICRS'], y=data_db['DE_ICRS'], color='purple', alpha=0.8)
plt.scatter(x=data_goudain['RA_ICRS'], y=data_goudain['DE_ICRS'], color='green', alpha=0.8)

print(type(data_db['DE_ICRS'][3]), type(data_db['RA_ICRS'][3]))
plt.show()
