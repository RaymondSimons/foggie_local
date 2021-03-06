import astropy
from astropy.io import fits
import numpy as np
from numpy import *
import math
from joblib import Parallel, delayed
import os, sys, argparse
import yt
from joblib import Parallel, delayed
from numpy import rec
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
    parser.add_argument('-simdir', '--simdir', default='/nobackupp2/mpeeples', help='simulation output directory')
    parser.add_argument('-figdir', '--figdir', default='/nobackupp2/rcsimons/foggie_momentum/figures/center_figures/central_gal', help='figures output directory')
    parser.add_argument('-haloname', '--haloname', default='halo_008508', help='halo_name')



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





def make_savefile(anchor_fits, simname,  haloname, simdir, DD, ds, ad):
    


    #a = fits.open(fits_name)


    '''
    mss = a['STAR_MASS'].data
    id_s = a['STARS_ID'].data
    xs, ys, zs = a['STARS_GAL_POSITION'].data
    vxs, vys, vzs = a['STARS_GAL_VELOCITY'].data
    ep_s = a['STARS_EPSILON_FIXED'].data
    xs_box, ys_box, zs_box = a['STARS_BOX_POSITION'].data
    vxs_box, vys_box, vzs_box = a['STARS_BOX_VELOCITY'].data
    '''
    DDname = 'DD%.4i'%DD
    #ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, simname,  DDname, DDname))
    #ad = ds.all_data()
    #def _stars(pfilter, data): return data[(pfilter.filtered_type, "particle_type")] == 2

    #yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
    #ds.add_particle_filter('stars')



    mss = ad['stars', 'particle_mass'].in_units('Msun')
    xs_box = ad['stars', 'particle_position_x'].in_units('kpc')
    ys_box = ad['stars', 'particle_position_y'].in_units('kpc')
    zs_box = ad['stars', 'particle_position_z'].in_units('kpc')  

    vxs_box = ad['stars', 'particle_velocity_x'].in_units('km/s')
    vys_box = ad['stars', 'particle_velocity_y'].in_units('km/s')
    vzs_box = ad['stars', 'particle_velocity_z'].in_units('km/s')

    id_s = ad['stars', 'particle_index']



    hdus = []
    prim_hdu = fits.PrimaryHDU()
    hdus.append(prim_hdu)

    for sat_n in arange(shape(anchor_fits)[0]):
        anchor_ids = anchor_fits[sat_n]
        gd_indices = []
        for anch_id in anchor_ids: 
            match = where(id_s == anch_id)[0]
            if len(match) > 0: gd_indices.append(int(match))
        gd_indices       = array(gd_indices)
 
        if len(gd_indices) > 10:
            print 'more than 10 anchor stars found for sat %i in DD%.4i'%(sat_n, DD)
            anchor_mss       = mss[gd_indices]
            anchor_xs_box    =  xs_box[gd_indices]
            anchor_ys_box    =  ys_box[gd_indices]
            anchor_zs_box    =  zs_box[gd_indices]
            anchor_vxs_box   = vxs_box[gd_indices]
            anchor_vys_box   = vys_box[gd_indices]
            anchor_vzs_box   = vzs_box[gd_indices]

            anchor_R = sqrt(anchor_xs_box**2. + anchor_ys_box**2. + anchor_zs_box**2.)

            hist_R, r_edges = np.histogram(anchor_R.value, weights = anchor_mss.value, bins = arange(min(anchor_R.value)-20,max(anchor_R.value)+20, 10))

            Rmid = np.mean([r_edges[argmax(hist_R)], r_edges[argmax(hist_R)+1]])
            good = where(abs(anchor_R.value - Rmid) < 5)[0]


            anchor_xs_box_avg, _  = weighted_avg_and_std(anchor_xs_box,  weights = anchor_mss, good = good)
            anchor_ys_box_avg, _  = weighted_avg_and_std(anchor_ys_box,  weights = anchor_mss, good = good)
            anchor_zs_box_avg, _  = weighted_avg_and_std(anchor_zs_box,  weights = anchor_mss, good = good)
            anchor_vxs_box_avg, _ = weighted_avg_and_std(anchor_vxs_box, weights = anchor_mss, good = good)
            anchor_vys_box_avg, _ = weighted_avg_and_std(anchor_vys_box, weights = anchor_mss, good = good)
            anchor_vzs_box_avg, _ = weighted_avg_and_std(anchor_vzs_box, weights = anchor_mss, good = good)

            box_avg = [anchor_xs_box_avg,
                       anchor_ys_box_avg,
                       anchor_zs_box_avg,
                       anchor_vxs_box_avg,
                       anchor_vys_box_avg,
                       anchor_vzs_box_avg]

            cols1 = fits.ColDefs([fits.Column(name = 'box_avg', array =  box_avg    , format = 'D'),
                                  fits.Column(name = 'anchor_mss     ', array =  anchor_mss    , format = 'D'),
                                  fits.Column(name = 'anchor_xs_box  ', array =  anchor_xs_box , format = 'D'),
                                  fits.Column(name = 'anchor_ys_box  ', array =  anchor_ys_box , format = 'D'),
                                  fits.Column(name = 'anchor_zs_box  ', array =  anchor_zs_box , format = 'D'),
                                  fits.Column(name = 'anchor_vxs_box ', array =  anchor_vxs_box, format = 'D'),
                                  fits.Column(name = 'anchor_vys_box ', array =  anchor_vys_box, format = 'D'),
                                  fits.Column(name = 'anchor_vzs_box ', array =  anchor_vzs_box, format = 'D'),
                                  fits.Column(name = 'ids_used_avg', array =  good, format = 'I'),
                                  ])
        else:
            print 'less than 10 anchor stars found for sat %i in DD%.4i'%(sat_n, DD)
            cols1 = fits.ColDefs([fits.Column(name = 'box_avg     ', array =  None   , format = '0D'),
                                  fits.Column(name = 'anchor_mss     ', array =  None, format = '0D'),
                                  fits.Column(name = 'anchor_xs_box  ', array =  None, format = '0D'),
                                  fits.Column(name = 'anchor_ys_box  ', array =  None, format = '0D'),
                                  fits.Column(name = 'anchor_zs_box  ', array =  None, format = '0D'),
                                  fits.Column(name = 'anchor_vxs_box ', array =  None, format = '0D'),
                                  fits.Column(name = 'anchor_vys_box ', array =  None, format = '0D'),
                                  fits.Column(name = 'anchor_vzs_box ', array =  None, format = '0D'),
                                  fits.Column(name = 'ids_used_avg', array =  None, format = '0D'),
                                  ])

        hdus.append(fits.BinTableHDU.from_columns(cols1, name = 'SAT_%.2i'%sat_n))

    hdus_fits = fits.HDUList(hdus)
    hdus_fits.writeto('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_anchorprops.fits'%(simname, DD), overwrite = True)




if __name__ == '__main__':
    args = parse()
    simname = args['simname']
    min_DD = int(args['DDmin'])
    max_DD = int(args['DDmax'])
    simdir = args['simdir']
    haloname = args['haloname']
    if 'selfshield_z15' in simname: anchor_simname = 'natural'
    else: anchor_simname = simname

    anchor_fits = np.load('/nobackupp2/rcsimons/foggie_momentum/catalogs/%s_anchors.npy'%anchor_simname)[()]

    def _stars(pfilter, data): return data[(pfilter.filtered_type, "particle_type")] == 2
    def run_tracker(DD, simdir, haloname, simname, anchor_fits):
        DDname = 'DD%.4i'%DD
        ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, simname,  DDname, DDname))
        ad = ds.all_data()

        yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
        ds.add_particle_filter('stars')
        make_savefile(anchor_fits = anchor_fits, simname = simname, haloname = haloname, simdir = simdir, DD = DD, ds = ds, ad = ad) 

    Parallel(n_jobs = 5)(delayed(run_tracker)(DD = DD, simdir = simdir, haloname = haloname, 
                                                                      simname = simname, anchor_fits = anchor_fits) for DD in np.arange(min_DD, max_DD))
    






















