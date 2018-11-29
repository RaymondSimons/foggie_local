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



if __name__ == '__main__':
    args = parse()

    simname = args['simname']
    DD = int(args['DD'])
    haloname = args['haloname']

    snapname = 'DD%.4i'%DD

    ds = yt.load('/nobackupp2/mpeeples/%s/%s/%s/%s'%(haloname, simname, snapname, snapname))

    dd = ds.all_data()
    star_mass = dd['stars', 'particle_mass'].in_units('Msun')

    n1 = len(star_mass)
    n2 = len(where(star_mass < 1.5e3)[0])
    z = ds.current_redshift


    np.save('/nobackupp2/rcsimons/foggie_momentum/nparticles/%s_DD%.4i_nparticles.npy'%(simname, DD), [n1, n2, z])





















