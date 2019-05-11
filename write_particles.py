#!/u/rcsimons/miniconda3/bin/python3.7
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





def make_savefile(simname,  haloname, simdir, DD, ds, ad):
    DDname = 'DD%.4i'%DD

    mss = ad['stars', 'particle_mass'].in_units('Msun')
    msd = ad['darkmatter', 'particle_mass'].in_units('Msun')

    ags =  ad['stars', 'age'].in_units('Gyr')
    agd =  ad['darkmatter', 'age'].in_units('Gyr')

    xs_box = ad['stars', 'particle_position_x'].in_units('kpc')
    ys_box = ad['stars', 'particle_position_y'].in_units('kpc')
    zs_box = ad['stars', 'particle_position_z'].in_units('kpc')  
    vxs_box = ad['stars', 'particle_velocity_x'].in_units('km/s')
    vys_box = ad['stars', 'particle_velocity_y'].in_units('km/s')
    vzs_box = ad['stars', 'particle_velocity_z'].in_units('km/s')

    xd_box = ad['darkmatter', 'particle_position_x'].in_units('kpc')
    yd_box = ad['darkmatter', 'particle_position_y'].in_units('kpc')
    zd_box = ad['darkmatter', 'particle_position_z'].in_units('kpc')  
    vxd_box = ad['darkmatter', 'particle_velocity_x'].in_units('km/s')
    vyd_box = ad['darkmatter', 'particle_velocity_y'].in_units('km/s')
    vzd_box = ad['darkmatter', 'particle_velocity_z'].in_units('km/s')

    ids = ad['stars', 'particle_index']
    idd = ad['darkmatter', 'particle_index']



    hdus = []
    prim_hdu = fits.PrimaryHDU()
    hdus.append(prim_hdu)

    colss = fits.ColDefs([fits.Column(name = 'id'       , array =  ids, format = 'D'),
                          fits.Column(name = 'mass'    ,  array =  mss    , format = 'D'),
                          fits.Column(name = 'x_box  ' ,  array =  xs_box , format = 'D'),
                          fits.Column(name = 'y_box  ' ,  array =  ys_box , format = 'D'),
                          fits.Column(name = 'z_box  ' ,  array =  zs_box , format = 'D'),
                          fits.Column(name = 'vx_box '  , array =  vxs_box, format = 'D'),
                          fits.Column(name = 'vy_box '  , array =  vys_box, format = 'D'),
                          fits.Column(name = 'vz_box '  , array =  vzs_box, format = 'D'),
                          fits.Column(name = 'age'      , array =  ags, format = 'D'),
                          ])

    colsd = fits.ColDefs([fits.Column(name = 'id'       , array =  idd, format = 'D'),
                          fits.Column(name = 'mass'    ,  array =  msd    , format = 'D'),
                          fits.Column(name = 'x_box  ' ,  array =  xd_box , format = 'D'),
                          fits.Column(name = 'y_box  ' ,  array =  yd_box , format = 'D'),
                          fits.Column(name = 'z_box  ' ,  array =  zd_box , format = 'D'),
                          fits.Column(name = 'vx_box '  , array =  vxd_box, format = 'D'),
                          fits.Column(name = 'vy_box '  , array =  vyd_box, format = 'D'),
                          fits.Column(name = 'vz_box '  , array =  vzd_box, format = 'D'),
                          fits.Column(name = 'age'      , array =  agd, format = 'D'),
                          ])



    hdus.append(fits.BinTableHDU.from_columns(colss, name = 'stars'))
    hdus.append(fits.BinTableHDU.from_columns(colsd, name = 'dark'))


    hdus_fits = fits.HDUList(hdus)
    hdus_fits.writeto('/nobackupp2/rcsimons/foggie_momentum/particles/%s/%s/%s_DD%.4i_particles.fits'%(haloname, simname, simname, DD), overwrite = True)




if __name__ == '__main__':
    args = parse()
    simname = args['simname']
    min_DD = int(args['DDmin'])
    max_DD = int(args['DDmax'])
    simdir = args['simdir']
    haloname = args['haloname']

    if simname == 'natural':      enzo_simname = 'natural'
    elif simname == 'natural_v2': enzo_simname = 'nref11n_v2_selfshield_z15'
    elif simname == 'natural_v3': enzo_simname = 'nref11n_v3_selfshield_z15'
    elif simname == 'natural_v4': enzo_simname = 'nref11n_v4_selfshield_z15'
    else: enzo_simname = simname

    def _stars(pfilter, data): return data[(pfilter.filtered_type, "particle_type")] == 2
    def _darkmatter(pfilter, data): return data[(pfilter.filtered_type, "particle_type")] == 4

    def run_writer(DD, simdir, haloname, simname):
        DDname = 'DD%.4i'%DD
        ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, enzo_simname,  DDname, DDname))
        ad = ds.all_data()

        yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
        yt.add_particle_filter("darkmatter",function=_darkmatter, filtered_type='all',requires=["particle_type"])
        ds.add_particle_filter('stars')
        ds.add_particle_filter('darkmatter')
        #if not os.path.exists('/nobackupp2/rcsimons/foggie_momentum/particles/%s/%s_DD%.4i_particles.fits'%(simname, simname, DD)):
        make_savefile(simname = simname, haloname = haloname, simdir = simdir, DD = DD, ds = ds, ad = ad) 

    Parallel(n_jobs = 1)(delayed(run_writer)(DD = DD, simdir = simdir, haloname = haloname, 
                                                                      simname = simname) for DD in np.arange(min_DD, max_DD))
    






















