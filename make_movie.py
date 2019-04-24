import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
import numpy as np
from numpy import *
import math
import seaborn as sns
from joblib import Parallel, delayed
import os, sys, argparse
import yt
yt.enable_parallelism()
from yt.units import kpc
from consistency import *
import seaborn as sns
from mpl_toolkits.axes_grid1 import AxesGrid

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




def make_figure(sat_n, figdir, wd, wdd, cen_fits, DD, cen_name, simdir, haloname, simname,  DDname):
       if len(cen_fits['SAT_%.2i'%sat_n].data['box_avg']) > 0:
            ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, simname,  DDname, DDname))
            def _stars(pfilter, data): return data[(pfilter.filtered_type, "particle_type")] == 2
            yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
            ds.add_particle_filter('stars')

            cenx = cen_fits['SAT_%.2i'%sat_n].data['box_avg'][0]
            ceny = cen_fits['SAT_%.2i'%sat_n].data['box_avg'][1]
            cenz = cen_fits['SAT_%.2i'%sat_n].data['box_avg'][2]
            cen_g = yt.YTArray([cenx, ceny, cenz], 'kpc')
            '''
            for axis in ['x', 'y', 'z']:

                if axis == 'x':
                    box = ds.r[cen_g[0] - 0.5 * yt.YTArray(wdd, 'kpc'): cen_g[0] + 0.5 * yt.YTArray(wdd, 'kpc'), \
                               cen_g[1] - 0.5 * yt.YTArray(5*wd,  'kpc'): cen_g[1] + 0.5 * yt.YTArray(5*wd,  'kpc'), \
                               cen_g[2] - 0.5 * yt.YTArray(5*wd,  'kpc'): cen_g[2] + 0.5 * yt.YTArray(wd,  'kpc')]

                elif axis == 'y':
                    box = ds.r[cen_g[0] - 0.5 * yt.YTArray(5*wd, 'kpc'): cen_g[0]   + 0.5 * yt.YTArray(5*wd, 'kpc'), \
                               cen_g[1] - 0.5 * yt.YTArray(wdd,  'kpc'): cen_g[1] + 0.5 * yt.YTArray(wdd,  'kpc'), \
                               cen_g[2] - 0.5 * yt.YTArray(5*wd,  'kpc'): cen_g[2]  + 0.5 * yt.YTArray(5*wd,  'kpc')]
                elif axis == 'z':
                    box = ds.r[cen_g[0] - 0.5 * yt.YTArray(5*wd, 'kpc'): cen_g[0]   + 0.5 * yt.YTArray(5*wd, 'kpc'), \
                               cen_g[1] - 0.5 * yt.YTArray(5*wd,  'kpc'): cen_g[1]  + 0.5 * yt.YTArray(5*wd,  'kpc'), \
                               cen_g[2] - 0.5 * yt.YTArray(wdd,  'kpc'): cen_g[2] + 0.5 * yt.YTArray(wdd,  'kpc')]




                fig = plt.figure(sat_n, figsize = (20,20))
                
                grid = AxesGrid(fig, (0.0,0.0,1.0,1.0),
                                nrows_ncols = (1, 2),
                                axes_pad = 0.0, label_mode = "1",
                                share_all = False, cbar_mode=None,
                                aspect = False)        



                W = yt.YTArray([wd, wd, wd], 'kpc')
                p = yt.ProjectionPlot(ds, axis, ("gas","density"), center = cen_g, data_source=box, width=W)
                p.set_unit(('gas','density'), 'Msun/pc**2')
                p.set_zlim(('gas', 'density'), zmin = density_proj_min * 0.1, zmax =  density_proj_max)
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
                
                figname = '%s_%.4i_%.2i_%s.png'%(cen_name, DD, sat_n, axis)
                fig.savefig('%s/%s'%(figdir,figname))
            '''

def do_plot(number):
    fig = plt.figure(number)

    a = random.sample(1000)
    b = random.sample(1000)

    # generate random data
    plt.scatter(a, b)

    plt.savefig("%03d.jpg" % (number,))
    plt.close()

    print("Done ", number)


if __name__ == '__main__':

    args = parse()
    simname = args['simname']
    DD = int(args['DD'])
    haloname = args['haloname']
    snapname = 'DD%.4i'%DD
    wd    = float(args['wd'])
    wdd    = float(args['wdd'])
    axis = args['axis']
    simdir = args['simdir']
    if simname == 'natural': cen_name = 'natural'
    if 'v2' in simname: cen_name = 'natural_v2'
    if 'v3' in simname: cen_name = 'natural_v3'
    if 'v4' in simname: cen_name = 'natural_v4'

    figdir = '/nobackupp2/rcsimons/foggie_momentum/sat_figures/%s'%cen_name
    if not os.path.exists(figdir): os.system('mkdir %s'%figdir)
    cen_fits = fits.open('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s/%s_DD%.4i_anchorprops.fits'%(cen_name, simname, DD))

    figs_list = []
    grid_list = []
    DDname = 'DD%.4i'%DD
    num_procs = 4

    #for sat_n in yt.parallel_objects(arange(3), num_procs):
    #    ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, simname,  DDname, DDname))
    #    def _stars(pfilter, data): return data[(pfilter.filtered_type, "particle_type")] == 2
    #    yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
    #    ds.add_particle_filter('stars')
    #    make_figure(sat_n, figdir, wd, wdd, ds, cen_fits, DD, cen_name)



    Parallel(n_jobs = 3)(delayed(make_figure)(sat_n, figdir, wd, wdd, cen_fits, DD, cen_name, simdir, haloname, simname, DDname) for sat_n in arange(3))


    '''
    for s, sat_n in enumerate(arange(1)):
        make_figure(sat_n, figdir, wd, wdd, ds, cen_fits, DD, cen_name)
        plt.close('all')
    '''
















