import healpy as hp
import healsparse as hsp
import numpy as np
from shapely.geometry.polygon import Polygon
import pandas as pd

# data file created by get_data.py
df = pd.read_csv('/global/u2/z/zdu863/notebooks/project/butler_all.csv',index_col=0)
# ra and dec of run 1.2i WFD field corners
wfd_ra = np.array([52.25, 52.11, 58.02, 57.87])
wfd_dec = np.array([-27.25, -32.25, -32.25, -27.25])

nsideCoverage = 32
nsideSparse = 2048

# get the the pixel indices within the DC2 1.2 field 
field_ang = np.array([np.radians(90-wfd_dec), np.radians(wfd_ra)])
ipix_infield = hp.query_polygon(nsideSparse, hp.ang2vec(*field_ang), nest=True, inclusive=True)
Nmpix = len(ipix_infield)

datanames = ['fiveSigmaDepth','airmass','rawSeeing','finSeeing','bg_mean','bg_var','zp','zp_err']
Nname = len(datanames)
bands = ['u','g','r','i','z','y']
nvisit = {b:np.zeros(Nmpix) for b in bands}
mp = {b:np.zeros([Nname, Nmpix]) for b in bands}

percent = 10
print('Computing systematic maps...')
# loop over all visits in the data file
for index, row in df.iterrows():
    if index/len(df)*100 > percent:
        print('%d%% done'%percent)
        percent += 10
    visit_ra = np.array([row['llcra'],row['lrcra'],row['urcra'],row['ulcra']]) #corner coordinates
    visit_dec = np.array([row['llcdec'],row['lrcdec'],row['urcdec'],row['ulcdec']])
    band = row['filter']
    frame_ang = np.array([np.radians(90-visit_dec), np.radians(visit_ra)]) #vertices of visit frame in theta_phi coords
    polyframe = Polygon(np.transpose(frame_ang)) #polygon of visit frame

    ipix_inframe = hp.query_polygon(nsideSparse, hp.ang2vec(*frame_ang), nest=True, inclusive=True) #pixel inside visit frame
    ipix_inter = set(ipix_infield).intersection(set(ipix_inframe)) #pixel inside visit frame and DC2 field

    if bool(ipix_inter):
        for ipix in ipix_inter:
            pixvec = hp.boundaries(nsideSparse,ipix,nest=True) #vertcies of healpix pixel
            polypix = Polygon(np.transpose(hp.vec2ang(np.transpose(pixvec)))) #polygon of healpix pixel
            if polyframe.contains(polypix):
                ind = np.where(ipix_infield==ipix)[0][0]
                nvisit[band][ind] += 1.
                for i,dataname in enumerate(datanames):
                    mp[band][i][ind] += row[dataname] * 1.
            elif polyframe.intersects(polypix):
                ind = np.where(ipix_infield==ipix)[0][0]
                area = polyframe.intersection(polypix).area/polypix.area #fraction of pixel inside visit frame
                nvisit[band][ind] += area
                for i,dataname in enumerate(datanames):
                    mp[band][i][ind] += row[dataname] * area

#write maps
dtype = [(dataname,'f8') for dataname in datanames]
dtype.insert(0,('nvisit','f8'))
for b in bands:
    norm_mp = np.zeros(Nmpix, dtype=dtype)
    norm_mp['nvisit'][nvisit[b]>0] = nvisit[b][nvisit[b]>0]
    for i,dataname in enumerate(datanames):
         norm_mp[dataname][nvisit[b]>0] = mp[b][i][nvisit[b]>0]/nvisit[b][nvisit[b]>0]
    hsp_mp = hsp.HealSparseMap.makeEmpty(nsideCoverage, nsideSparse, dtype=dtype, primary='nvisit')
    hsp_mp.updateValues(ipix_infield, norm_mp)
    hsp_mp.write('/global/u2/z/zdu863/notebooks/project/fitsmap/allmaps-' + b + '.fits', clobber=True)