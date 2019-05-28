import GCRCatalogs
import numpy as np
import healsparse as hsp
import healpy as hp
from GCR import GCRQuery

def depth_map(ra, dec, mags, nsideSparse, field_corners):
    good = np.logical_or(np.logical_not(np.isnan(ra)),np.logical_not(np.isnan(dec)))
    pix_infield = hp.query_polygon(nsideSparse, hp.ang2vec(*field_corners), nest=True, inclusive=True)
    pix_nums = hp.ang2pix(nsideSparse, np.radians(90-dec), np.radians(ra), nest=True) # pixels where the objects are located

    map_out = np.zeros(len(pix_infield)) + hp.UNSEEN
    map_var_out = np.zeros(len(pix_infield)) + hp.UNSEEN
    
    for pix in np.unique(pix_nums):
        if pix in pix_infield:
            ind = np.where(pix_infield==pix)[0][0]
        else:
            continue
        mask = pix==pix_nums # mask of all objects inside the pixel
        if len(mags[mask])>0:
            map_out[ind] = np.mean(mags[mask])
            map_var_out[ind] = np.std(mags[mask])
    
    dtype = [('out','float'), ('var_out','float')]
    rec_mp = np.rec.fromarrays([map_out, map_var_out], dtype=dtype)
    hsp_mp = hsp.HealSparseMap.makeEmpty(32, nsideSparse, dtype=dtype, primary='out')
    hsp_mp.updateValues(pix_infield, rec_mp)
    
    return hsp_mp

catalog = GCRCatalogs.load_catalog('dc2_object_run1.2i_all_columns_with_photoz')
band = 'i'
snr_cuts = [
    GCRQuery('clean'), 
    GCRQuery('detect_isPrimary'),
    GCRQuery('snr_%s_cModel > 4'%band), # SNR > 4
    GCRQuery('snr_%s_cModel < 6'%band), # SNR < 6
    GCRQuery((np.isfinite, 'mag_%s_cModel'%band)),
    GCRQuery((np.isfinite, 'cModelFluxErr_%s'%band)),
    GCRQuery((np.isfinite, 'cModelFlux_%s'%band))
]

# Loads the data after cut
data_cut = catalog.get_quantities(['ra', 'dec', 'mag_%s_cModel'%band], 
                              filters = snr_cuts)
ra,dec = data_cut['ra'], data_cut['dec']
mags = data_cut['mag_%s_cModel'%band]

nsideSparse = 2048
wfd_ra = np.array([52.25, 52.11, 58.02, 57.87])
wfd_dec = np.array([-27.25, -32.25, -32.25, -27.25])
field_ang = np.array([np.radians(90-wfd_dec), np.radians(wfd_ra)])
hsp_mp = depth_map(ra, dec, mags, nsideSparse, field_ang)
hsp_mp.write('depth_map.fits', clobber=True)
