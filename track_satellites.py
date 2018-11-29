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



def weighted_avg_and_std(values, weights, good):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    values = values[good]
    weights = weights[good]
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))





def make_savefile(anchor_fits, DD, simname):
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



    hdus = []
    prim_hdu = fits.PrimaryHDU()
    hdus.append(prim_hdu)

    for sat_n in arange(5):

        anchor_ids = anchor_fits['OLDSTARS_%.2i'%sat_n].data['ids']
        gd_indices = array([0 for i in arange(len(anchor_ids))])




        for g in arange(len(anchor_ids)):
            gd_indices[g] = int(where(id_s == anchor_ids[g])[0])

        anchor_mss      = mss[gd_indices]
        anchor_xs       = xs[gd_indices]
        anchor_ys       = ys[gd_indices]
        anchor_zs       = zs[gd_indices]
        anchor_vxs      = vxs[gd_indices]
        anchor_vys      = vys[gd_indices]
        anchor_vzs      = vzs[gd_indices]

        anchor_xs_box    =  xs_box[gd_indices]
        anchor_ys_box    =  ys_box[gd_indices]
        anchor_zs_box    =  zs_box[gd_indices]
        anchor_vxs_box   = vxs_box[gd_indices]
        anchor_vys_box   = vys_box[gd_indices]
        anchor_vzs_box   = vzs_box[gd_indices]


        good = where((abs(anchor_xs - median(anchor_xs)) < 5) & (abs(anchor_ys - median(anchor_ys)) < 5) & (abs(anchor_zs - median(anchor_zs)) < 5))          





        anchor_xs_avg, _ = weighted_avg_and_std(anchor_xs,  weights = anchor_mss, good = good)
        anchor_ys_avg, _ = weighted_avg_and_std(anchor_ys,  weights = anchor_mss, good = good)
        anchor_zs_avg, _ = weighted_avg_and_std(anchor_zs,  weights = anchor_mss, good = good)
        anchor_vxs_avg, _= weighted_avg_and_std(anchor_vxs, weights = anchor_mss, good = good)
        anchor_vys_avg, _= weighted_avg_and_std(anchor_vys, weights = anchor_mss, good = good)
        anchor_vzs_avg, _= weighted_avg_and_std(anchor_vzs, weights = anchor_mss, good = good)
        anchor_xs_box_avg, _  = weighted_avg_and_std(anchor_xs_box,  weights = anchor_mss, good = good)
        anchor_ys_box_avg, _  = weighted_avg_and_std(anchor_ys_box,  weights = anchor_mss, good = good)
        anchor_zs_box_avg, _  = weighted_avg_and_std(anchor_zs_box,  weights = anchor_mss, good = good)
        anchor_vxs_box_avg, _ = weighted_avg_and_std(anchor_vxs_box, weights = anchor_mss, good = good)
        anchor_vys_box_avg, _ = weighted_avg_and_std(anchor_vys_box, weights = anchor_mss, good = good)
        anchor_vzs_box_avg, _ = weighted_avg_and_std(anchor_vzs_box, weights = anchor_mss, good = good)

        cols1 = fits.Column(name = 'anchor_mss     ', array =  anchor_mss    , format = 'D')
         fits.Column(name = 'anchor_xs      ', array =  anchor_xs     , format = 'D')
         fits.Column(name = 'anchor_ys      ', array =  anchor_ys     , format = 'D')
         fits.Column(name = 'anchor_zs      ', array =  anchor_zs     , format = 'D')
         fits.Column(name = 'anchor_vxs     ', array =  anchor_vxs    , format = 'D')
         fits.Column(name = 'anchor_vys     ', array =  anchor_vys    , format = 'D')
         fits.Column(name = 'anchor_vzs     ', array =  anchor_vzs    , format = 'D')
         fits.Column(name = 'anchor_xs_box  ', array =  anchor_xs_box , format = 'D')
         fits.Column(name = 'anchor_ys_box  ', array =  anchor_ys_box , format = 'D')
         fits.Column(name = 'anchor_zs_box  ', array =  anchor_zs_box , format = 'D')
         fits.Column(name = 'anchor_vxs_box ', array =  anchor_vxs_box, format = 'D')
         fits.Column(name = 'anchor_vys_box ', array =  anchor_vys_box, format = 'D')
         fits.Column(name = 'anchor_vzs_box ', array =  anchor_vzs_box, format = 'D')



        hdus.append(fits.BinTableHDU.from_columns(cols1, name = 'table_%.2i'%sat_n))





        to_save = [anchor_xs_box_avg, anchor_ys_box_avg, anchor_zs_box_avg, anchor_vxs_box_avg, anchor_vys_box_avg, anchor_vzs_box_avg, anchor_xs_avg, anchor_ys_avg, anchor_zs_avg, anchor_vxs_avg, anchor_vys_avg, anchor_vzs_avg]



        print '/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_sat%.2i_cen.npy'%(simname, DD, sat_n)
        np.save('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_sat%.2i_cen.npy'%(simname, DD, sat_n),to_save)


    hdus_fits = fits.HDUList(hdus)
    hdus_fits.writeto('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_sat%.2i_anchorprops.fits'%(simname, DD, sat_n), overwrite = True)




if __name__ == '__main__':
    args = parse()

    simname = args['simname']
    min_DD = int(args['DDmin'])
    max_DD = int(args['DDmax'])
    momentum_directory = '/nobackupp2/rcsimons/foggie_momentum/momentum_fits'    
    #anchor_ids = np.load('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_anchors.npy'%simname)

    anchor_fits = fits.open('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_anchors_DD0250.fits'%simname)

    Parallel(n_jobs = -1, backend = 'threading')(delayed(make_savefile)(anchor_fits = anchor_fits, DD = DD, simname = simname) for DD in np.arange(min_DD, max_DD))

























