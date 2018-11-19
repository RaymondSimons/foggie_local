import sys
import os
import glob
import yt
import numpy as np
from numpy import *
import astropy
from astropy.cosmology import Planck13 as cosmo
import findGalaxyProps as fGP
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

    parser.add_argument('-run_parallel', '--run_parallel', default=False, help='Run parallel')

    parser.add_argument('-simname', '--simname', default=None, help='Simulation to be analyzed.')

    parser.add_argument('-snapname', '--snapname', default=None, help='Snapshot files to be analyzed.')

    parser.add_argument('-haloname', '--haloname', default='halo_008508', help='halo_name')
    parser.add_argument('-ddmin', '--ddmin', default=500, help='halo_name')
    parser.add_argument('-ddmax', '--ddmax', default=600, help='halo_name')
    parser.add_argument('-n_jobs', '--n_jobs', default=3, help='number of jobs')



    args = vars(parser.parse_args())
    return args





def get_particles(haloname, simname, snapname):
    snaps = np.sort(np.asarray(glob.glob("/nobackupp2/mpeeples/%s/%s/%s/%s"%(haloname, simname, snapname, snapname))))


    #abssnap = os.path.abspath(snaps[0])
    assert os.path.lexists(snaps[0])

    out_dir = '/u/rcsimons/'

    assert os.path.lexists(out_dir)

    new_snapfiles = np.asarray(snaps)

    ts = yt.DatasetSeries(new_snapfiles)




    def _stars(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 2

    def _darkmatter(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 4

    yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
    yt.add_particle_filter("darkmatter",function=_darkmatter, filtered_type='all',requires=["particle_type"])

    for ds,snapfile in zip(reversed(ts),np.flipud(new_snapfiles)):
        ad = ds.all_data()
        ds.add_particle_filter('stars')
        ds.add_particle_filter('darkmatter')

        dark_pos_x  = ad['darkmatter', 'particle_position_x'].in_units('kpc')
        stars_pos_x = ad['stars', 'particle_position_x'].in_units('kpc')

        print shape(dark_pos_x), shape(stars_pos_x)



if __name__ == "__main__":

    args = parse()

    simname = args['simname']
    snapname = args['snapname']
    haloname = args['haloname']
    run_parallel = args['run_parallel']

    print haloname, simname, snapname, run_parallel
    snaps = np.sort(np.asarray(glob.glob("/nobackupp2/mpeeples/%s/%s/%s/%s"%(haloname, simname, snapname, snapname))))


    #abssnap = os.path.abspath(snaps[0])
    assert os.path.lexists(snaps[0])

    out_dir = '/u/rcsimons/'

    assert os.path.lexists(out_dir)

    new_snapfiles = np.asarray(snaps)

    ts = yt.DatasetSeries(new_snapfiles)




    def _stars(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 2

    #this gets dark matter particles in zoom region only
    def _darkmatter(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 4

    yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
    yt.add_particle_filter("darkmatter",function=_darkmatter, filtered_type='all',requires=["particle_type"])

    for ds,snapfile in zip(reversed(ts),np.flipud(new_snapfiles)):
        ad = ds.all_data()
        ds.add_particle_filter('stars')
        ds.add_particle_filter('darkmatter')

        dark_pos_x  = ad['darkmatter', 'particle_position_x'].in_units('kpc')
        stars_pos_x = ad['stars', 'particle_position_x'].in_units('kpc')
        all_particles= dd['io', 'particle_position_x'].in_units('kpc')
        print shape(dark_pos_x), shape(stars_pos_x)













