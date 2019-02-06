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
from consistency import *
import matplotlib.pyplot as plt
import yt
import numpy as np
from yt.units import kpc
import matplotlib.pyplot as plt
from consistency import *
import seaborn as sns
from mpl_toolkits.axes_grid1 import AxesGrid
plt.ioff()

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



def make_figure(sat_n, figname, figdir, wd, wdd, ds):

        try:
            cen_np = np.load('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_sat%.2i_1049_cen.npy'%(simname, DD, sat_n))[()]
        except:
            return
        cenx = cen_np[0]
        ceny = cen_np[1]
        cenz = cen_np[2]
        cen_g = yt.YTArray([cenx, ceny, cenz], 'kpc')



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





            fig = plt.figure(1, figsize = (20,20))

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
            fig.savefig('%s/satn%i_%s-axis_%s'%(figdir, sat_n, axis, figname))
            plt.close('all')

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
    figdir = args['figdir']
    figname = args['figname']



    figs_list = []
    grid_list = []




    DDname = 'DD%.4i'%DD
    ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, simname,  DDname, DDname))

    def _stars(pfilter, data): return data[(pfilter.filtered_type, "particle_type")] == 2

    yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
    ds.add_particle_filter('stars')

    #Parallel(n_jobs = 5, backend = 'threading')(delayed(make_figure)(sat_n, figname, figdir, wd, wdd, ds, figs_list[s], grid_list[s]) for s, sat_n in enumerate(arange(5)))



    for s, sat_n in enumerate(arange(5)):

        make_figure(sat_n, figname, figdir, wd, wdd, ds)
        plt.close('all')




































'''





#Ls = [Lx, Ly, Lz]
Ls = [0, 0, 1]
ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, simname, snapname, snapname))

if True:
    cenx = float(args['cenx'])
    ceny = float(args['ceny'])
    cenz = float(args['cenz'])    
else:
    #Read center from galprops
    galaxy_props_file = galprops_outdir + '/'  + simname + '_' + snapname + '_galprops.npy'
    gp = np.load(galaxy_props_file)[()]
    cenx = float(gp['stars_center'][0])
    ceny = float(gp['stars_center'][1])
    cenz = float(gp['stars_center'][2])







cen = yt.YTArray([cenx, ceny, cenz], 'kpc')

density_proj_min = 5e-2  # msun / pc^2
density_proj_max = 1e4

metal_min = 1.e-4
metal_max = 10.



print Ls

#p = yt.off_axis_projection(ds, cen, Ls, W, N, ('gas', 'density'), north_vector =  north_vector)#, zmin = density_proj_min, zmax = density_proj_max)
#M = yt.off_axis_projection(ds, cen, Ls, W, N, ('gas', 'metallicity'), weight = ('gas', 'density'), north_vector =  north_vector)#, zmin = density_proj_min, zmax = density_proj_max)


p = yt.ProjectionPlot(ds, fields = 'density', center = cen, width =  W)#, zmin = density_proj_min, zmax = density_proj_max)
M = yt.ProjectionPlot(ds, fields = 'metallicity', center = cen, width =  W, weight = 'density', north_vector =  north_vector)#, zmin = density_proj_min, zmax = density_proj_max)

density_color_map = sns.blend_palette(("black", "#4575b4", "#4daf4a", "#ffe34d", "darkorange"), as_cmap=True)
metal_color_map = sns.blend_palette(("black", "#4575b4", "#984ea3", "#984ea3", "#d73027", "darkorange", "#ffe34d"), as_cmap=True)

p = p.in_units('Msun * pc**-2')
M = M.in_units('Zsun')


if not add_cbar:
    fig, axes = plt.subplots(1,2, figsize = (20,10))
    fig.subplots_adjust(left = 0.0, right = 1.0, top =1.0, bottom = 0.0, hspace = 0.0, wspace = 0.0)

else:
    fig, axes = plt.subplots(1,1, figsize = (10,9))
    fig.subplots_adjust(left = 0.0, right = 0.92, top =1.0, bottom = 0.0, hspace = 0.0, wspace = 0.0)



im1 = axes[0].imshow(np.log10(p), vmin = log10(density_proj_min), vmax = log10(density_proj_max), cmap = density_color_map)
im2 = axes[1].imshow(np.log10(M), vmin = log10(metal_min),        vmax = log10(metal_max), cmap = metal_color_map)

for ax in axes:
    ax.axis('off')

if add_cbar:
    cax = fig.add_axes([0.915, 0.0, 0.02, 1.0])
    cbr = fig.colorbar(im1, cax=cax, orientation="vertical", cmap = density_color_map)
    cbr.set_label(r'log projected gas density (M$_{\odot}$ kpc$^{-2}$)', fontsize = 15)

fig.savefig("%s/%s_%.4i.png"%(figdir, simname, DD), dpi = 300)

'''
















