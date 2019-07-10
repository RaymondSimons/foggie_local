#!/u/rcsimons/miniconda3/bin/python3.7
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy import *
import numpy as np
from numpy import *
import math
import seaborn as sns
from joblib import Parallel, delayed
import os, sys, argparse
import yt
from yt.units import kpc
from consistency import *
import seaborn as sns
from mpl_toolkits.axes_grid1 import AxesGrid
import warnings
import copy
warnings.filterwarnings("ignore")
import time
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
    parser.add_argument('-DDmin', '--DDmin', default=600, help='min DD to use')
    parser.add_argument('-DDmax', '--DDmax', default=600, help='max DD to use')

    parser.add_argument('-simname', '--simname', default=None, help='Simulation to be analyzed.')

    parser.add_argument('-haloname', '--haloname', default='halo_008508', help='halo_name')

    parser.add_argument('-cenx', '--cenx', default=None, help='box position of galaxy, x')
    parser.add_argument('-ceny', '--ceny', default=None, help='box position of galaxy, y')
    parser.add_argument('-cenz', '--cenz', default=None, help='box position of galaxy, z')
    parser.add_argument('-wd', '--wd', default=yt.YTArray([75, 75, 75], 'kpc'), help='width of camera, kpc')
    parser.add_argument('-wdd', '--wdd', default=100., help='width of camera, kpc')
    parser.add_argument('-simdir', '--simdir', default='/nobackupp2/mpeeples', help='simulation output directory')
    parser.add_argument('-figdir', '--figdir', \
                        default='/nobackupp2/rcsimons/foggie_momentum/figures/satellites', \
                        help='figures output directory')
    parser.add_argument('-figname', '--figname', default='temp.png', help='figures output directory')
    parser.add_argument('-axis', '--axis', default='z', help='figures output directory')


    args = vars(parser.parse_args())
    return args




def make_figure(figdir, DD, cen_name, simdir, haloname, simname,  wd = 50., wdd = 100., wd2 = 150., wdd2 = 150., wd3 = 1000., wdd3 = 1000.):
        DDname = 'DD%.4i'%DD
        ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, simname,  DDname, DDname))

        cen_fits = np.load('/nobackupp2/rcsimons/foggie_momentum/catalogs/sat_interpolations/%s_interpolations_DD0150_new.npy'%cen_name, allow_pickle=True)[()]



        cenx = cen_fits['CENTRAL']['fxe'](DD)
        ceny = cen_fits['CENTRAL']['fye'](DD)
        cenz = cen_fits['CENTRAL']['fze'](DD)

        cen_g = yt.YTArray([cenx, ceny, cenz], 'kpc')



        box = ds.r[cen_g[0] - 0.5 * yt.YTArray(3*wd, 'kpc'): cen_g[0]   + 0.5 * yt.YTArray(3*wd, 'kpc'), \
                   cen_g[1] - 0.5 * yt.YTArray(3*wd,  'kpc'): cen_g[1]  + 0.5 * yt.YTArray(3*wd,  'kpc'), \
                   cen_g[2] - 0.5 * yt.YTArray(wdd,  'kpc'): cen_g[2] + 0.5 * yt.YTArray(wdd,  'kpc')]



        axis = 'z'

        density_proj_min = 2.5e-1
        p = yt.ProjectionPlot(ds, axis, ("gas","density"), center = cen_g, data_source=box, width=(W, 'kpc'))
        p.set_unit(('gas','density'), 'Msun/pc**2')
        p.set_zlim(('gas', 'density'), zmin = density_proj_min, zmax =  density_proj_max)
        p.set_cmap(('gas', 'density'), density_color_map)
        p.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
        p.hide_axes()
        p.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
        p.annotate_scale(size_bar_args={'color':'white'})
        p.save('%s/%.4i_density.png'%(figdir, DD), mpl_kwargs = {'dpi': 500})

        metal_color_map = sns.blend_palette(
            ("black", "#4575b4", "#5d31c4", "#984ea3", "#d73027",
             "darkorange", "#ffe34d"), as_cmap=True)


        metal_min = 2.e-3
        p = yt.ProjectionPlot(ds, axis, "metallicity", center = cen_g, data_source=box, width=(W, 'kpc'), weight_field = ('gas', 'density'))
        p.set_unit(('gas','metallicity'), 'Zsun')
        p.set_zlim(('gas', 'metallicity'), zmin = metal_min, zmax =  metal_max)
        p.set_cmap(('gas', 'metallicity'), metal_color_map)
        p.set_log("metallicity", True)
        p.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
        p.hide_axes()
        p.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
        p.annotate_scale(size_bar_args={'color':'white'})
        p.save('%s/%.4i_metallicity.png'%(figdir, DD), mpl_kwargs = {'dpi': 500})



if __name__ == '__main__':

    args = parse()
    simname = args['simname']
    DDmin = int(args['DDmin'])
    DDmax = int(args['DDmax'])
    haloname = args['haloname']
    axis = args['axis']
    simdir = args['simdir']
    if simname == 'natural': cen_name = 'natural'
    if 'v2' in simname: cen_name = 'natural_v2'
    if 'v3' in simname: cen_name = 'natural_v3'
    if 'v4' in simname: cen_name = 'natural_v4'
    if simname == 'nref11n_nref10f': cen_name = 'nref11n_nref10f'    
    if simname == 'nref11c_nref9f': cen_name = 'nref11c_nref9f'    

    figdir = '/nobackupp2/rcsimons/foggie_momentum/figures/density_metallicity_movies'
    if not os.path.exists(figdir): os.system('mkdir %s'%figdir)

    lst = []
    for DD in arange(DDmin, DDmax):
        make_figure(figdir, DD, cen_name, simdir, haloname, simname)        

