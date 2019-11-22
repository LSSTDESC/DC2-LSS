#!/usr/bin/env python

#download data for star-galaxy classification
import sys
sys.path.insert(0, '/global/common/software/lsst/common/miniconda/py3-4.5.12/envs/stack/lib/python3.7/site-packages')

import os
os.chdir('/global/cscratch1/sd/zdu863/objcat_sgclass')

import GCRCatalogs
from GCR import GCRQuery
import re
import glob
import numpy as np
import pandas as pd

# flags 
standard_cuts =[
    GCRQuery('detect_isPrimary'),
    GCRQuery('modelfit_CModel_flag_badCentroid==False'),
    GCRQuery('base_SdssCentroid_flag==False'),
    GCRQuery('deblend_skipped==False'),
    GCRQuery('base_PixelFlags_flag_edge==False'),
    GCRQuery('base_PixelFlags_flag_interpolatedCenter==False'),
    GCRQuery('base_PixelFlags_flag_saturatedCenter==False'),
    GCRQuery('base_PixelFlags_flag_crCenter==False'),
    GCRQuery('base_PixelFlags_flag_bad==False'),
    GCRQuery('base_PixelFlags_flag_suspectCenter==False'),
    GCRQuery('base_NaiveCentroid_flag==False')
]

# additional cuts under reduce_cat.py#L151
blend_cut = [
    GCRQuery('base_Blendedness_abs<=0.42169650342')
]

# SNR>3 in ugizy-bands and SNR>6 in r-band 
ugizyflux_cut = GCRQuery('+'.join('(snr_{}_cModel > 3)*1'.format(b) for b in 'ugizy') + '>= 2')
rflux_cut = GCRQuery('snr_r_cModel > 6')
flux_cut = [ugizyflux_cut & rflux_cut]

tracts = np.arange(int(sys.argv[1]),int(sys.argv[2])) # range of tract indexes where we download the data from
path = "./sg_tract_*.csv"
temp = [2897] + [re.findall(r'\d+', filename)[0] for filename in glob.glob(path)]
existed_tracts = np.array(temp).astype(int) # tracts that's already downloaded
remained_tracts = tracts[~np.isin(tracts, existed_tracts)] # tracts still need to be downloaded

if len(remained_tracts)==0:
    print(tracts[0],'to',tracts[-1],'is done.')
else:
    remained_tracts = np.append(remained_tracts, np.amax(existed_tracts)) # starts downloading from the largest existing tract
    catalog = GCRCatalogs.load_catalog('dc2_object_run2.1i_dr1b_all_columns')
    for tract in remained_tracts:
        data = catalog.get_quantities(['ra', 'dec', 'mag_r_cModel', 'extendedness'],
                                      filters = standard_cuts + blend_cut + flux_cut + 
                                      [GCRQuery((np.isfinite, 'mag_r_cModel')), 
                                       GCRQuery((np.isfinite, 'extendedness'))],
                                      native_filters=['tract==%d'%tract])
        
        df = pd.DataFrame(data)
        df.to_csv('./sg_tract_%d.csv'%tract)
        