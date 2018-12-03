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

    parser.add_argument('-cenx', '--cenx', default=None, help='box position, x')
    parser.add_argument('-ceny', '--ceny', default=None, help='box position, y')
    parser.add_argument('-cenz', '--cenz', default=None, help='box position, z')
    parser.add_argument('-ls', '--ls', default=[[1,0,0], [0,1,0], [0, 0, 1]], help='direction of camera')
    parser.add_argument('-w', '--w', default=yt.YTArray([100, 100, 100], 'kpc'), help='width of camera, kpc')
    parser.add_argument('-n', '--n', default=[0,0.7,0.7], help='north vector of camera')
    parser.add_argument('-npix', '--npix', default=512, help='number of pixels')

    args = vars(parser.parse_args())
    return args



if __name__ == '__main__':
    args = parse()
    simname = args['simname']
    DD = int(args['DD'])
    haloname = args['haloname']
    snapname = 'DD%.4i'%DD
    cenx = float(args['cenx'])
    ceny = float(args['ceny'])
    cenz = float(args['cenz'])
    Ls   = args['ls']
    W    = args['w']
    north_vector = args['n']
    N = int(args['npix'])
    Ls = [[1,0,0], [0,1,0], [0, 0, 1]]

    W1 = yt.YTArray([15, 15, 15], 'kpc')
    W2 = yt.YTArray([100, 100, 100], 'kpc')

    north_vector = [0,0.7,0.7]
    N = 512

    ds = yt.load('/nobackupp2/mpeeples/%s/%s/%s/%s'%(haloname, simname, snapname, snapname))
    cen = yt.YTArray([cenx, ceny, cenz], 'kpc')


    image1 = yt.off_axis_projection(ds, cen, Ls, W, N, ('gas', 'density'), north_vector =  north_vector)
    yt.write_image(np.log10(image1), "/nobackupp2/rcsimons/foggie_momentum/figures/center_figures/%s_offaxis_projection.png"%simname)





#python make_movie_frames.py -DD 300 -simname nref11n_selfshield_z15 -cenx 10000 -ceny 10000 -cenz 10000 
















