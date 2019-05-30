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




def make_figure(figdir, DD, cen_name, simdir, haloname, simname,  wd = 100., wdd = 100., wd2 = 150., wdd2 = 150., wd3 = 1000., wdd3 = 1000.):
        #figname_zoomoutfar = '%s_%.4i_%.2i_%s_zoomoutfar.png'%(cen_name, DD, 6, 'x')
        #figname_check = '%s/%s/%s/%s'%(figdir,'x', 'zoomoutfar', figname_zoomoutfar)
        #if os.path.isfile(figname_check): return

        DDname = 'DD%.4i'%DD
        ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, simname,  DDname, DDname))
        def _stars(pfilter, data): return data[(pfilter.filtered_type, "particle_type")] == 2
        yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
        ds.add_particle_filter('stars')
        cen_fits = np.load('/nobackupp2/rcsimons/foggie_momentum/catalogs/sat_interpolations/%s_interpolations_DD0150_new.npy'%cen_name, allow_pickle=True)[()]

        '''
        #encoding='latin1' is needed for loading python 2 pickles in python 3
        if 'natural' in cen_name:
            central_xyz_fit = np.load('/nobackupp2/rcsimons/foggie_momentum/catalogs/center_natural.npy', allow_pickle=True, encoding='latin1')[()]

        elif ('nref11n_nref10f' in cen_name) | ('nref11c_nref9f' in cen_name):
            central_xyz_fit = np.load('/nobackupp2/rcsimons/foggie_momentum/catalogs/center_nref11n_nref10f.npy', allow_pickle=True, encoding='latin1')[()]
        '''

        #central_x = cen_fits['CENTRAL']['fxe'](DD)
        #central_y = cen_fits['CENTRAL']['fye'](DD)
        #central_z = cen_fits['CENTRAL']['fze'](DD)

        central_x = cen_fits['CENTRAL']['fxe'](DD)
        central_y = cen_fits['CENTRAL']['fye'](DD)
        central_z = cen_fits['CENTRAL']['fze'](DD)

        print (central_x, central_y, central_z)

        '''
        xf = central_xyz_fit['x']
        yf = central_xyz_fit['y']
        zf = central_xyz_fit['z']
        '''
        #central_x = xf[0] * DD**4. + xf[1] * DD**3. + xf[2] * DD**2. + xf[3] * DD + xf[4]
        #central_y = yf[0] * DD**4. + yf[1] * DD**3. + yf[2] * DD**2. + yf[3] * DD + yf[4]
        #central_z = zf[0] * DD**4. + zf[1] * DD**3. + zf[2] * DD**2. + zf[3] * DD + zf[4]

        cen_central = yt.YTArray([central_x, central_y, central_z], 'kpc')


        W = yt.YTArray([wd, wd, wd], 'kpc')
        W2 = yt.YTArray([wd2, wd2, wd2], 'kpc')
        W3 = yt.YTArray([wd3, wd3, wd3], 'kpc')


        for axis in ['x']:#@, 'y', 'z']:
            if axis == 'x':                
                box2 = ds.r[cen_central[0] - 0.5 * yt.YTArray(wdd2, 'kpc'):   cen_central[0] + 0.5 * yt.YTArray(wdd2, 'kpc'), \
                            cen_central[1] - 0.5 * yt.YTArray(3*wd2,  'kpc'): cen_central[1] + 0.5 * yt.YTArray(3*wd2,  'kpc'), \
                            cen_central[2] - 0.5 * yt.YTArray(3*wd2,  'kpc'): cen_central[2] + 0.5 * yt.YTArray(3*wd2,  'kpc')]


                box3 = ds.r[cen_central[0] - 0.5 * yt.YTArray(wdd3, 'kpc'):   cen_central[0] + 0.5 * yt.YTArray(wdd3, 'kpc'), \
                            cen_central[1] - 0.5 * yt.YTArray(3*wd3,  'kpc'): cen_central[1] + 0.5 * yt.YTArray(3*wd3,  'kpc'), \
                            cen_central[2] - 0.5 * yt.YTArray(3*wd3,  'kpc'): cen_central[2] + 0.5 * yt.YTArray(3*wd3,  'kpc')]


            if axis == 'y':

                box2 = ds.r[cen_central[0] - 0.5 * yt.YTArray(3*wd2, 'kpc') : cen_central[0]  + 0.5 * yt.YTArray(3*wd2, 'kpc'), \
                            cen_central[1] - 0.5 * yt.YTArray(wdd2,  'kpc') : cen_central[1]  + 0.5 * yt.YTArray(wdd2,  'kpc'), \
                            cen_central[2] - 0.5 * yt.YTArray(3*wd2,  'kpc'): cen_central[2]  + 0.5 * yt.YTArray(3*wd2,  'kpc')]


                box3 = ds.r[cen_central[0] - 0.5 * yt.YTArray(3*wd3, 'kpc') : cen_central[0]  + 0.5 * yt.YTArray(3*wd3, 'kpc'), \
                            cen_central[1] - 0.5 * yt.YTArray(wdd3,  'kpc') : cen_central[1]  + 0.5 * yt.YTArray(wdd3,  'kpc'), \
                            cen_central[2] - 0.5 * yt.YTArray(3*wd3,  'kpc'): cen_central[2]  + 0.5 * yt.YTArray(3*wd3,  'kpc')]


            if axis == 'z':

                box2 = ds.r[cen_central[0] - 0.5 * yt.YTArray(3*wd2, 'kpc'):  cen_central[0]   + 0.5 * yt.YTArray(3*wd2, 'kpc'), \
                            cen_central[1] - 0.5 * yt.YTArray(3*wd2,  'kpc'): cen_central[1]  + 0.5 * yt.YTArray(3*wd2,  'kpc'), \
                            cen_central[2] - 0.5 * yt.YTArray(wdd2,  'kpc'):  cen_central[2] + 0.5 * yt.YTArray(wdd2,  'kpc')]


                box3 = ds.r[cen_central[0] - 0.5 * yt.YTArray(3*wd3, 'kpc'):  cen_central[0]   + 0.5 * yt.YTArray(3*wd3, 'kpc'), \
                            cen_central[1] - 0.5 * yt.YTArray(3*wd3,  'kpc'): cen_central[1]  + 0.5 * yt.YTArray(3*wd3,  'kpc'), \
                            cen_central[2] - 0.5 * yt.YTArray(wdd3,  'kpc'):  cen_central[2] + 0.5 * yt.YTArray(wdd3,  'kpc')]

            '''

            a = time.time()
            p_wd2_g= yt.ProjectionPlot(ds, axis, ("gas","density"), center = cen_central, data_source=box2, width=W2)
            b = time.time()
            print ('p_wd2_g', b-a)

            p_wd2_g.set_unit(('gas','density'), 'Msun/pc**2')
            p_wd2_g.set_zlim(('gas', 'density'), zmin = density_proj_min, zmax =  density_proj_max)
            p_wd2_g.set_cmap(('gas', 'density'), density_color_map)
            p_wd2_g.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
            p_wd2_g.hide_axes()
            p_wd2_g.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
            p_wd2_g.annotate_scale(size_bar_args={'color':'white'})



            a = time.time()
            p_wd2_s = yt.ParticleProjectionPlot(ds, axis, ('stars', 'particle_mass'), center = cen_central, data_source=box2, width = W2)   
            b = time.time()
            print ('p_wd2_s', b-a)

            cmp = plt.cm.Greys_r
            cmp.set_bad('k')
            p_wd2_s.set_cmap(field = ('stars','particle_mass'), cmap = cmp)
            p_wd2_s.hide_axes()
            p_wd2_s.annotate_scale(size_bar_args={'color':'white'})

            p_wd2_s.set_zlim(field = ('stars','particle_mass'), zmin = 2.e35 * 0.3, zmax = 1.e42*0.9)



            a = time.time()
            p_wd3_g= yt.ProjectionPlot(ds, axis, ("gas","density"), center = cen_central, data_source=box3, width=W3)
            b = time.time()
            print ('p_wd3_g', b-a)

            p_wd3_g.set_unit(('gas','density'), 'Msun/pc**2')
            p_wd3_g.set_zlim(('gas', 'density'), zmin = density_proj_min, zmax =  density_proj_max)
            p_wd3_g.set_cmap(('gas', 'density'), density_color_map)
            p_wd3_g.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
            p_wd3_g.hide_axes()
            p_wd3_g.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
            p_wd3_g.annotate_scale(size_bar_args={'color':'white'})


            a = time.time()
            p_wd3_s = yt.ParticleProjectionPlot(ds, axis, ('stars', 'particle_mass'), center = cen_central, data_source=box3, width = W3)   
            b = time.time()
            print ('p_wd3_s', b-a)
        
    
            cmp = plt.cm.Greys_r
            cmp.set_bad('k')
            p_wd3_s.set_cmap(field = ('stars','particle_mass'), cmap = cmp)
            p_wd3_s.hide_axes()
            p_wd3_s.annotate_scale(size_bar_args={'color':'white'})
            p_wd3_s.set_zlim(field = ('stars','particle_mass'), zmin = 2.e35 * 0.3, zmax = 1.e42*0.9)
            '''

            for sat_n in arange(6,7):
                if sat_n < 6:
                    cenx = cen_fits['SAT_%.2i'%sat_n]['fxe'](DD)
                    ceny = cen_fits['SAT_%.2i'%sat_n]['fye'](DD)
                    cenz = cen_fits['SAT_%.2i'%sat_n]['fze'](DD)

                if sat_n == 6:
                    #cenx = cen_fits['CENTRAL']['fxe'](DD)
                    #ceny = cen_fits['CENTRAL']['fye'](DD)
                    #cenz = cen_fits['CENTRAL']['fze'](DD)

                    cenx = cen_fits['CENTRAL']['fxe'](DD)
                    ceny = cen_fits['CENTRAL']['fye'](DD)
                    cenz = cen_fits['CENTRAL']['fze'](DD)


                cen_g = yt.YTArray([cenx, ceny, cenz], 'kpc')
                print ('satellite center: ', cen_g)
                print ('central center: ', cen_central)

                figname_zoomin  = '%s_%.4i_%.2i_%s_zoomin_100kpc.png'%(cen_name, DD, sat_n, axis)
                figname_zoomout = '%s_%.4i_%.2i_%s_zoomout.png'%(cen_name, DD, sat_n, axis)
                figname_zoomoutfar = '%s_%.4i_%.2i_%s_zoomoutfar.png'%(cen_name, DD, sat_n, axis)

                if axis == 'x':
                    box = ds.r[cen_g[0] - 0.5 * yt.YTArray(wdd, 'kpc'): cen_g[0]   + 0.5 * yt.YTArray(wdd, 'kpc'), \
                               cen_g[1] - 0.5 * yt.YTArray(3*wd,  'kpc'): cen_g[1] + 0.5 * yt.YTArray(3*wd,  'kpc'), \
                               cen_g[2] - 0.5 * yt.YTArray(3*wd,  'kpc'): cen_g[2] + 0.5 * yt.YTArray(3*wd,  'kpc')]

                    p_1 = cen_g[1] - cen_central[1] 
                    p_2 = cen_g[2] - cen_central[2]
                    p_3 = cen_g[0] - cen_central[0]

                elif axis == 'y':
                    box = ds.r[cen_g[0] - 0.5 * yt.YTArray(3*wd, 'kpc'): cen_g[0]   + 0.5 * yt.YTArray(3*wd, 'kpc'), \
                               cen_g[1] - 0.5 * yt.YTArray(wdd,  'kpc'): cen_g[1]   + 0.5 * yt.YTArray(wdd,  'kpc'), \
                               cen_g[2] - 0.5 * yt.YTArray(3*wd,  'kpc'): cen_g[2]  + 0.5 * yt.YTArray(3*wd,  'kpc')]
                    p_1 = cen_g[2] - cen_central[2] 
                    p_2 = cen_g[0] - cen_central[0]
                    p_3 = cen_g[1] - cen_central[1]

                elif axis == 'z':
                    box = ds.r[cen_g[0] - 0.5 * yt.YTArray(3*wd, 'kpc'): cen_g[0]   + 0.5 * yt.YTArray(3*wd, 'kpc'), \
                               cen_g[1] - 0.5 * yt.YTArray(3*wd,  'kpc'): cen_g[1]  + 0.5 * yt.YTArray(3*wd,  'kpc'), \
                               cen_g[2] - 0.5 * yt.YTArray(wdd,  'kpc'): cen_g[2] + 0.5 * yt.YTArray(wdd,  'kpc')]

                    p_1 = cen_g[0] - cen_central[0] 
                    p_2 = cen_g[1] - cen_central[1]
                    p_3 = cen_g[2] - cen_central[2]


                fig = plt.figure(sat_n)
                
                grid = AxesGrid(fig, (0.0,0.0,1.0,1.0),
                                nrows_ncols = (1, 2),
                                axes_pad = 0.0, label_mode = "1",
                                share_all = False, cbar_mode=None,
                                aspect = False)        

                p = yt.ProjectionPlot(ds, axis, ("gas","density"), center = cen_g, data_source=box, width=W)
                p.set_unit(('gas','density'), 'Msun/pc**2')
                p.set_zlim(('gas', 'density'), zmin = density_proj_min, zmax =  density_proj_max)
                p.set_cmap(('gas', 'density'), density_color_map)
                p.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
                p.hide_axes()
                p.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
                p.annotate_scale(size_bar_args={'color':'white'})
                plot = p.plots[("gas","density")]
                plot.figure = fig
                plot.axes = grid[0].axes
                p._setup_plots()


                p = yt.ParticleProjectionPlot(ds, axis, ('stars', 'particle_mass'), center = cen_g, data_source=box, width = W)   
                cmp = plt.cm.Greys_r
                cmp.set_bad('k')
                p.set_cmap(field = ('stars','particle_mass'), cmap = cmp)
                p.hide_axes()
                p.annotate_scale(size_bar_args={'color':'white'})

                p.set_zlim(field = ('stars','particle_mass'), zmin = 2.e35 * 0.3, zmax = 1.e42*0.9)
                plot = p.plots[('stars','particle_mass')]
                plot.figure = fig
                plot.axes = grid[1].axes
                p._setup_plots()

                fig.set_size_inches(12, 6)
                fig.savefig('%s/%s/%s/%s'%(figdir,axis, 'zoomin', figname_zoomin))
                plt.close(fig)

                '''
                fig = plt.figure(sat_n)
                
                grid = AxesGrid(fig, (0.0,0.0,1.0,1.0),
                                nrows_ncols = (1, 2),
                                axes_pad = 0.0, label_mode = "1",
                                share_all = False, cbar_mode=None,
                                aspect = False)        

                p = copy.copy(p_wd2_g)
                plot = p.plots[("gas","density")]
                plot.figure = fig
                plot.axes = grid[0].axes
                p._setup_plots()
                print (abs(p_1), abs(p_2), W2)
                if (abs(p_1) < wd2/2.) & (abs(p_2) < wd2/2.) & (abs(p_3) < wd2/2.): plot.axes.scatter(p_1, p_2, marker = 'o', facecolor = "none", edgecolor='red', lw = 2, s = 800)

                p = copy.copy(p_wd2_s)
                plot = p.plots[('stars','particle_mass')]
                plot.figure = fig
                plot.axes = grid[1].axes
                p._setup_plots()
                #if (abs(p_1) < wd2/2.) & (abs(p_2) < wd2/2.) & (abs(p_3) < wd2/2.): plot.axes.scatter(p_1, p_2, marker = 'o', facecolor = "none", edgecolor='red', lw = 2, s = 800)
                fig.set_size_inches(12, 6)                    
                fig.savefig('%s/%s/%s/%s'%(figdir,axis, 'zoomout', figname_zoomout))
                plt.close(fig)





                fig = plt.figure(sat_n)
                
                grid = AxesGrid(fig, (0.0,0.0,1.0,1.0),
                                nrows_ncols = (1, 2),
                                axes_pad = 0.0, label_mode = "1",
                                share_all = False, cbar_mode=None,
                                aspect = False)        
                p = copy.copy(p_wd3_g)
                plot = p.plots[("gas","density")]
                plot.figure = fig
                plot.axes = grid[0].axes
                p._setup_plots()
                print (abs(p_1), abs(p_2), W2)
                #if (abs(p_1) < wd3/2.) & (abs(p_2) < wd3/2.) & (abs(p_3) < wdd3/2.): plot.axes.scatter(p_1, p_2, marker = 'o', facecolor = "none", edgecolor='red', lw = 2, s = 800)

                p = copy.copy(p_wd3_s)
                plot = p.plots[('stars','particle_mass')]
                plot.figure = fig
                plot.axes = grid[1].axes
                p._setup_plots()
                #if (abs(p_1) < wd3/2.) & (abs(p_2) < wd3/2.) & (abs(p_3) < wdd3/2.): plot.axes.scatter(p_1, p_2, marker = 'o', facecolor = "none", edgecolor='red', lw = 2, s = 800)
                fig.set_size_inches(12, 6)                    
                fig.savefig('%s/%s/%s/%s'%(figdir,axis, 'zoomoutfar', figname_zoomoutfar))
                plt.close(fig)
                '''

        ds.index.clear_all_data()        

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

    figdir = '/nobackupp2/rcsimons/foggie_momentum/sat_figures/%s'%cen_name
    if not os.path.exists(figdir): os.system('mkdir %s'%figdir)

    lst = []
    for DD in arange(DDmin, DDmax):
        make_figure(figdir, DD, cen_name, simdir, haloname, simname)        

