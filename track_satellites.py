import astropy
from astropy.io import fits
import numpy as np
from numpy import *
import math
from joblib import Parallel, delayed
import os, sys, argparse


def parse():
    '''
    Parse command line arguments
    ''' 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''\
                                Generate the cameras to use in Sunrise and make projection plots
                                of the data for some of these cameras. Then export the data within
                                the fov to a FITS file in a format that Sunrise understands.
                                ''')
    parser.add_argument('-DDmax', '--DDmax', default=None, help='max DD')
    parser.add_argument('-DDmin', '--DDmin', default=None, help='min DD')
    parser.add_argument('-simname', '--simname', default=None, help='Simulation to be analyzed.')



    args = vars(parser.parse_args())
    return args



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

    print 'Opening %s...'%fits_name

    a = fits.open(fits_name)



    mss = a['STAR_MASS'].data
    id_s = a['STARS_ID'].data
    xs, ys, zs = a['STARS_GAL_POSITION'].data
    vxs, vys, vzs = a['STARS_GAL_VELOCITY'].data
    ep_s = a['STARS_EPSILON_FIXED'].data
    xs_box, ys_box, zs_box = a['STARS_BOX_POSITION'].data
    vxs_box, vys_box, vzs_box = a['STARS_BOX_VELOCITY'].data




    for i in arange(5):

        anchor_ids = a['OLDSTARS_%.2i'%i].data['ids']
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

        anchor_xs_avg, _ = weighted_avg_and_std(anchor_xs,  weights = anchor_mss)
        anchor_ys_avg, _ = weighted_avg_and_std(anchor_ys,  weights = anchor_mss)
        anchor_zs_avg, _ = weighted_avg_and_std(anchor_zs, weights = anchor_mss)
        anchor_vxs_avg, _= weighted_avg_and_std(anchor_vxs, weights = anchor_mss)
        anchor_vys_avg, _= weighted_avg_and_std(anchor_vys, weights = anchor_mss)
        anchor_vzs_avg, _= weighted_avg_and_std(anchor_vzs, weights = anchor_mss)
        anchor_xs_box_avg, _  = weighted_avg_and_std(anchor_xs_box, weights = anchor_mss)
        anchor_ys_box_avg, _  = weighted_avg_and_std(anchor_ys_box, weights = anchor_mss)
        anchor_zs_box_avg, _  = weighted_avg_and_std(anchor_zs_box, weights = anchor_mss)
        anchor_vxs_box_avg, _ = weighted_avg_and_std(anchor_vxs_box, weights = anchor_mss)
        anchor_vys_box_avg, _ = weighted_avg_and_std(anchor_vys_box, weights = anchor_mss)
        anchor_vzs_box_avg, _ = weighted_avg_and_std(anchor_vzs_box, weights = anchor_mss)


        to_save = [anchor_xs_box_avg, anchor_ys_box_avg, anchor_zs_box_avg, anchor_vxs_box_avg, anchor_vys_box_avg, anchor_vzs_box_avg, anchor_xs_avg, anchor_ys_avg, anchor_zs_avg, anchor_vxs_avg, anchor_vys_avg, anchor_vzs_avg]



        print '/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_sat%.2i_cen.npy'%(simname, DD, i)
        np.save('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_sat%.2i_cen.npy'%(simname, DD, i),to_save)




if __name__ == '__main__':
    args = parse()

    simname = args['simname']
    min_DD = int(args['DDmin'])
    max_DD = int(args['DDmax'])
    momentum_directory = '/nobackupp2/rcsimons/foggie_momentum/momentum_fits'    
    #anchor_ids = np.load('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_anchors.npy'%simname)

    anchor_fits = fits.open('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_anchors_DD0250.fits'%simname)

    Parallel(n_jobs = -1, backend = 'threading')(delayed(make_savefile)(anchor_ids = anchor_fits, DD = DD, simname = simname) for DD in np.arange(min_DD, max_DD))

























