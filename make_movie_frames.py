import matplotlib
matplotlib.use('Agg')
import astropy
from astropy.io import fits
import numpy as np
from numpy import *
import math
import seaborn as sns
from joblib import Parallel, delayed
import os, sys, argparse
import yt
import matplotlib.pyplot as plt
plt.close('all')
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
    parser.add_argument('-DD', '--DD', default=600, help='DD to use')

    parser.add_argument('-simname', '--simname', default=None, help='Simulation to be analyzed.')

    parser.add_argument('-haloname', '--haloname', default='halo_008508', help='halo_name')

    parser.add_argument('-cenx', '--cenx', default=None, help='box position, x')
    parser.add_argument('-ceny', '--ceny', default=None, help='box position, y')
    parser.add_argument('-cenz', '--cenz', default=None, help='box position, z')
    parser.add_argument('-ls', '--ls', default=[1,0,0], help='direction of camera')
    parser.add_argument('-w', '--w', default=yt.YTArray([85, 85, 85], 'kpc'), help='width of camera, kpc')
    parser.add_argument('-n', '--n', default=[0,0.7,0.7], help='north vector of camera')
    parser.add_argument('-npix', '--npix', default=512, help='number of pixels')
    parser.add_argument('-simdir', '--simdir', default='/nobackupp2/mpeeples', help='simulation output directory')
    parser.add_argument('-figdir', '--figdir', default='/nobackupp2/rcsimons/foggie_momentum/figures/center_figures', help='figures output directory')
    parser.add_argument('-add_cbar', '--add_cbar', default=False, help='add a colorbar')

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
    simdir = args['simdir']
    figdir = args['figdir']
    add_cbar = args['add_cbar']


    ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, simname, snapname, snapname))
    cen = yt.YTArray([cenx, ceny, cenz], 'kpc')

    density_proj_min = 5e-2  # msun / pc^2
    density_proj_max = 1e4

    p = yt.off_axis_projection(ds, cen, Ls, W, N, ('gas', 'density'), north_vector =  north_vector)#, zmin = density_proj_min, zmax = density_proj_max)
    density_color_map = sns.blend_palette(("black", "#4575b4", "#4daf4a", "#ffe34d", "darkorange"), as_cmap=True)

    p = p.in_units('Msun * pc**-2')


    if not add_cbar:
        fig, ax = plt.subplots(1,1, figsize = (10,10))
        fig.subplots_adjust(left = 0.0, right = 1.0, top =1.0, bottom = 0.0)

    else:
        fig, ax = plt.subplots(1,1, figsize = (10,9))
        fig.subplots_adjust(left = 0.0, right = 0.92, top =1.0, bottom = 0.0, hspace = 0.0, wspace = 0.0)



    im1 = ax.imshow(np.log10(p), vmin = log10(density_proj_min), vmax = log10(density_proj_max), cmap = density_color_map)
    ax.axis('off')
    if add_cbar:
        cax = fig.add_axes([0.915, 0.0, 0.02, 1.0])
        cbr = fig.colorbar(im1, cax=cax, orientation="vertical", cmap = density_color_map)
        cbr.set_label(r'log projected gas density (M$_{\odot}$ kpc$^{-2}$)', fontsize = 15)

    fig.savefig("%s/%s_%.4i.png"%(figdir, simname, DD), dpi = 300)


















