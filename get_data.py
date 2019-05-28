import warnings
warnings.filterwarnings("ignore")
import sqlite3
from sqlite3 import Error
import numpy as np
import pandas as pd
import lsst.daf.persistence as dp
import astropy.wcs

data_imsim = '/global/cscratch1/sd/desc/DC2/data/Run1.2i_globus_in2p3_20181217/w_2018_39/rerun/281118' # Path to the processed images
butler = dp.Butler(data_imsim)
datarefs = butler.subset('src').cache # This contains information about the visits that have been processed
corners = np.array([[-0.5,-0.5],[-0.5,3999.5],[4071.5,-0.5],[4071.5,3999.5]])

# query using butler
query = {'visit':[],'filter':[],'bg_mean':[],'bg_var':[],'zp':[],'zp_err':[],'llcra':[],'llcdec':[],'lrcra':[],'lrcdec':[],'ulcra':[],'ulcdec':[]
        ,'urcra':[],'urcdec':[]}
# loop over all visits
for data in datarefs:
    try:
        metadata = butler.get('calexp_md', data).toDict()
    except Exception:
        continue
    query['visit'].append(data['visit'])
    query['filter'].append(data['filter'])
    query['bg_mean'].append(metadata['BGMEAN'])
    query['bg_var'].append(metadata['BGVAR'])
    query['zp'].append(metadata['FLUXMAG0'])
    query['zp_err'].append(metadata['FLUXMAG0ERR'])
    w = astropy.wcs.WCS(metadata)
    c = w.wcs_pix2world(corners,0)
    query['llcra'].append(c[0,0])
    query['llcdec'].append(c[0,1])
    query['ulcra'].append(c[1,0])
    query['ulcdec'].append(c[1,1])
    query['lrcra'].append(c[2,0])
    query['lrcdec'].append(c[2,1])
    query['urcra'].append(c[3,0])
    query['urcdec'].append(c[3,1])    
df = pd.DataFrame(query)

def create_connection(db_file):
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
    return None

def get_all_visits(conn):
    cur = conn.cursor()
    cur.execute("SELECT obsHistID, airmass, altitude, ditheredRA, ditheredDec, filtSkyBrightness, finSeeing, rawSeeing, visitExpTime, fiveSigmaDepth FROM ObsHistory")
    rows = cur.fetchall()
    return rows

# query the database
db_file = '/global/projecta/projectdirs/lsst/groups/SSim/DC2/minion_1016_desc_dithered_v4.db'
conn = create_connection(db_file)
query = {'visit':[],'airmass':[],'altitude':[],'ra':[],'dec':[],'filtSkyBrightness':[],'finSeeing':[],'rawSeeing':[],'visitExpTime':[],'fiveSigmaDepth':[]}
data = get_all_visits(conn)
convert = {'altitude','ra','dec'}
for item in data:
    for key,value in zip(query.keys(),item):
        if key in convert:
            value = np.rad2deg(value)
        query[key].append(value)        
df2 = pd.DataFrame(query)

# Create a new dataframe to store the list of visits that appear in the v1.2i processed images (access by butler)
# This is because the visits in the processed images are a subset of those in the database, which contains 
# conditions for the full survey
newdict = {key:[] for key in df2.columns.values[1:]}
for i in range(len(df)):
    visit = df.iloc[i]['visit']
    for key in newdict.keys():
        newdict[key].append(df2.loc[df2['visit']==visit][key].values[0])
df3 = pd.DataFrame(newdict)

#combine data from butler and database queries
df4 = df.join(df3)
df4.to_csv('butler_all.csv')