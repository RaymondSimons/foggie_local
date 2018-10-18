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









if __name__ == "__main__":

    args = parse()

    simname = args['simname']
    snapname = args['snapname']
    haloname = args['haloname']
    run_parallel = args['run_parallel']

    print haloname, simname, snapname, run_parallel
