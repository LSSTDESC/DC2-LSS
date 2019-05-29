import sys
sys.path.insert(0, '/global/common/software/lsst/common/miniconda/current/lib/python3.6/site-packages')

import GCRCatalogs
import numpy as np
import healsparse as hsp
import healpy as hp
from GCR import GCRQuery

def depth_map_meanSNRrange(ra, dec, mags, snr, snrthreshold, nsideSparse):
    # 5sigma depth= mean of mags of galaxies with 4<SNR<6
    pix_nums = hp.ang2pix(nsideSparse, np.radians(90-dec), np.radians(ra), nest=True)
    pix_uni =np.unique(pix_nums)
    
    map_out = np.zeros(len(pix_uni)) + hp.UNSEEN
    map_var_out = np.zeros(len(pix_uni)) + hp.UNSEEN
    map_nobj = np.zeros(len(pix_uni)) + hp.UNSEEN
    
    # For each healpix pixel    
    for ind, pix in enumerate(pix_uni):
        # Select all objects within this pixel
        mask = pix==pix_nums
        # Compute mean of magnitudes with in snr restriction
        maskedmags = mags[mask][(snr[mask]>snrthreshold-1) & (snr[mask]<snrthreshold+1)]
        if len(maskedmags)>0:
            map_out[ind] = np.mean(maskedmags)
            map_var_out[ind] = np.std(maskedmags)
    
    dtype = [('out','float'), ('var_out','float')]
    rec_mp = np.rec.fromarrays([map_out, map_var_out], dtype=dtype)
    hsp_mp = hsp.HealSparseMap.makeEmpty(32, nsideSparse, dtype=dtype, primary='out')
    hsp_mp.updateValues(pix_uni, rec_mp)
    
    return hsp_mp

catalog = GCRCatalogs.load_catalog('dc2_object_run1.2i_all_columns_with_photoz')
band = 'i'
simple_cuts = [
    GCRQuery('clean'), 
    GCRQuery('detect_isPrimary'),
    GCRQuery((np.isfinite, 'ra')),
    GCRQuery((np.isfinite, 'dec')),
    GCRQuery((np.isfinite, 'mag_%s_cModel'%band)),
    GCRQuery((np.isfinite, 'snr_%s_cModel'%band))
]

# Loads the data after cut
data_cut = catalog.get_quantities(['ra', 'dec', 'snr_%s_cModel'%band, 'mag_%s_cModel'%band], 
                                  filters = simple_cuts)
ra,dec = data_cut['ra'], data_cut['dec']
mags = data_cut['mag_%s_cModel'%band]
snr = data_cut['snr_%s_cModel'%band]
hsp_mp = depth_map_meanSNRrange(ra, dec, mags, snr, 5, 2048)
hsp_mp.write('depth_map.fits', clobber=True)