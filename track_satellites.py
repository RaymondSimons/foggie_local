import astropy
from astropy.io import fits
import numpy as np
from numpy import *
import math
from joblib import Parallel, delayed


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))





def make_savefile(anchor_ids, DD, simname):
    fits_name = momentum_directory + '/' + simname + '_' + 'DD%.4i_momentum.fits'%DD
    a = fits.open(fits_name)
    mss = a['STAR_MASS'].data
    id_s = a['STARS_ID'].data
    xs, ys, zs = a['STARS_GAL_POSITION'].data
    vxs, vys, vzs = a['STARS_GAL_VELOCITY'].data
    ep_s = a['STARS_EPSILON_FIXED'].data
    xs_box, ys_box, zs_box = a['STARS_BOX_POSITION'].data
    vxs_box, vys_box, vzs_box = a['STARS_BOX_VELOCITY'].data


    gd_indices = array([0 for i in arange(len(anchor_ids))])
    for g in arange(len(anchor_ids)):
        gd_indices[g] = int(where(id_s == anchor_ids[g])[0])

    anchor_mss    = mss[gd_indices]
    anchor_xs    = xs[gd_indices]
    anchor_ys    = ys[gd_indices]
    anchor_zs    = zs[gd_indices]
    anchor_vxs   = vxs[gd_indices]
    anchor_vys   = vys[gd_indices]
    anchor_vzs   = vzs[gd_indices]

    anchor_xs_box    =  xs_box[gd_indices]
    anchor_ys_box    =  ys_box[gd_indices]
    anchor_zs_box    =  zs_box[gd_indices]
    anchor_vxs_box   = vxs_box[gd_indices]
    anchor_vys_box   = vys_box[gd_indices]
    anchor_vzs_box   = vzs_box[gd_indices]


    anchor_xs_box_avg, _  = weighted_avg_and_std(anchor_xs_box, weights = anchor_mss)
    anchor_ys_box_avg, _  = weighted_avg_and_std(anchor_ys_box, weights = anchor_mss)
    anchor_zs_box_avg, _  = weighted_avg_and_std(anchor_zs_box, weights = anchor_mss)
    anchor_vxs_box_avg, _ = weighted_avg_and_std(anchor_vxs_box, weights = anchor_mss)
    anchor_vys_box_avg, _ = weighted_avg_and_std(anchor_vys_box, weights = anchor_mss)
    anchor_vzs_box_avg, _ = weighted_avg_and_std(anchor_vzs_box, weights = anchor_mss)


    to_save = [anchor_xs_box_avg, anchor_ys_box_avg, anchor_zs_box_avg, anchor_vxs_box_avg, anchor_vys_box_avg, anchor_vzs_box_avg]

    print 'Saving to /nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_cen.npy...'%(simname, DD)
    np.save('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_cen.npy'%(simname, DD),to_save)




if __name__ == '__main__':
    min_DD = 600
    max_DD = 610
    momentum_directory = '/nobackupp2/rcsimons/foggie_momentum/momentum_fits'    
    for simname in ['nref11n_selfshield_z15']:
        anchor_ids = np.load('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_anchors.npy'%simname)
        Parallel(n_jobs = 3, backend = 'threading')(delayed(make_savefile)(anchor_ids = anchor_ids, DD = DD, simname = simname) for DD in np.arange(min_DD, max_DD))

























