import matplotlib
matplotlib.use('Agg')
import astropy
from astropy.io import fits
import numpy as np
from numpy import *
import math
from joblib import Parallel, delayed
import os, sys, argparse
import yt
import matplotlib.pyplot as plt

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
    parser.add_argument('-DD', '--DD', default=None, help='DD to use')

    parser.add_argument('-simname', '--simname', default=None, help='Simulation to be analyzed.')

    parser.add_argument('-haloname', '--haloname', default='halo_008508', help='halo_name')

    args = vars(parser.parse_args())
    return args



def measure_mass(simname, DD, sat_n):
    cen_file =  np.load('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_sat%.2i_cen.npy'%(simname, DD, sat_n))[()]



    anchor_xs_box_avg, anchor_ys_box_avg, anchor_zs_box_avg, anchor_vxs_box_avg, anchor_vys_box_avg, anchor_vzs_box_avg, anchor_xs_avg, anchor_ys_avg, anchor_zs_avg, anchor_vxs_avg, anchor_vys_avg, anchor_vzs_avg = cen_file

    cen = yt.YTArray([anchor_xs_box_avg, anchor_ys_box_avg, anchor_zs_box_avg], 'kpc')


    gc_sphere =  ds.sphere(cen, ds.arr(5,'kpc'))



    DM_mass = gc_sphere.quantities.total_quantity([("darkmatter", "particle_mass")]).to('Msun')
    gas_mass = gc_sphere.quantities.total_quantity([("gas", "cell_mass")]).to('Msun')
    gas_metal_mass = gc_sphere.quantities.total_quantity([("gas", "metal_mass")]).to('Msun')
    stars_mass = gc_sphere.quantities.total_quantity([("stars", "particle_mass")]).to('Msun')
    youngstars_mass = gc_sphere.quantities.total_quantity([("youngstars", "particle_mass")]).to('Msun')


    mass = [gas_mass,gas_metal_mass, DM_mass, stars_mass, youngstars_mass]

    np.save('/nobackupp2/rcsimons/foggie_momentum/satellite_masses/%s_DD%.4i_mass_sat%.2i.npy'%(simname, DD, sat_n), mass)





if __name__ == '__main__':
    args = parse()
    simname = args['simname']
    DD = int(args['DD'])
    haloname = args['haloname']
    snapname = 'DD%.4i'%DD
    ds = yt.load('/nobackupp2/mpeeples/%s/%s/%s/%s'%(haloname, simname, snapname, snapname))


    def _stars(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 2

    def _youngstars(pfilter, data):
        return data[(pfilter.filtered_type, "age")] < 2.e7

    # these are only the must refine dark matter particles
    def _darkmatter(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 4

    yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
    yt.add_particle_filter("youngstars",function=_youngstars, filtered_type='all',requires=["age"])
    yt.add_particle_filter("darkmatter",function=_darkmatter, filtered_type='all',requires=["particle_type"])

    ds.add_particle_filter('stars')
    ds.add_particle_filter('darkmatter')
    ds.add_particle_filter('youngstars')


    Parallel(n_jobs = -1, backend = 'threading')(delayed(measure_mass)(simname = simname, DD = DD, sat_n = sat_n, ds = ds) for sat_n in np.arange(5))






























