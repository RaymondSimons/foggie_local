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
    parser.add_argument('-lx', '--lx', default=1, help='direction of camera, x')
    parser.add_argument('-ly', '--ly', default=0, help='direction of camera, y')
    parser.add_argument('-lz', '--lz', default=0, help='direction of camera, z')
    parser.add_argument('-w', '--w', default=yt.YTArray([50, 50, 50], 'kpc'), help='width of camera, kpc')
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
    Lx   = float(args['lx'])
    Ly   = float(args['ly'])
    Lz   = float(args['lz'])
    cenx   = float(args['cenx'])
    ceny   = float(args['ceny'])
    cenz   = float(args['cenz'])
    W    = args['w']
    north_vector = args['n']
    N = int(args['npix'])
    simdir = args['simdir']
    figdir = args['figdir']


    wd = 150
    W = yt.YTArray([wd, wd, wd], 'kpc')


    ds = yt.load('%s/%s/%s/%s'%(simdir, simname, DDname, DDname))
    gp = np.load('/Users/rsimons/Dropbox/rcs_foggie/galprops/halo_008508/nref11n_nref10f_selfshield_z6_DD1049_galprops.npy', encoding = 'latin1')[()]

    cen_g = yt.YTArray([cenx, ceny, cenz], 'kpc')



    sats_n = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/nref11n_selfshield_z15_%s_sats.npy'%DDname)[()]
    sats_f = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/nref11n_nref10f_selfshield_z6_%s_sats.npy'%DDname)[()]
    def _stars(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 2

    # these are only the must refine dark matter particles
    def _darkmatter(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 4

    yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
    yt.add_particle_filter("darkmatter",function=_darkmatter, filtered_type='all',requires=["particle_type"])
    ds_n.add_particle_filter('stars')
    ds_n.add_particle_filter('darkmatter')
    ds_f.add_particle_filter('stars')
    ds_f.add_particle_filter('darkmatter')




    fig = plt.figure(1, figsize = (20,20))

    grid = AxesGrid(fig, (0.0,0.0,1.0,1.0),
                    nrows_ncols = (1, 2),
                    axes_pad = 0.0, label_mode = "1",
                    share_all = False, cbar_mode=None,
                    aspect = False)        

    box = ds.r[cen_g[0] - 0.5 * yt.YTArray(max([200., wd]), 'kpc'): cen_g[0] + 0.5 * yt.YTArray(max([200., wd]), 'kpc'), \
               cen_g[1] - 0.5 * yt.YTArray(max([200., wd]), 'kpc'): cen_g[1] + 0.5 * yt.YTArray(max([200., wd]), 'kpc'), \
               cen_g[2] - 0.5 * yt.YTArray(max([200., wd]), 'kpc'): cen_g[2] + 0.5 * yt.YTArray(max([200., wd]), 'kpc')]



    p = yt.ProjectionPlot(ds, 'y', ("gas","density"), center = cen_g, data_source=box, width=W)
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


    p = yt.ParticleProjectionPlot(ds, 'y', ('stars', 'particle_mass'), center = cen_g, data_source=box, width = W)   
    cmp = plt.cm.Greys_r
    cmp.set_bad('k')
    p.set_cmap(field = ('stars','particle_mass'), cmap = cmp)
    p.hide_axes()
    p.annotate_scale(size_bar_args={'color':'white'})
    #p.annotate_marker(pos = ds.arr(cen), coord_system = 'data', marker = 'o')
    #c = [0.1, 0.1]
    #p.annotate_text(c, 'o', coord_system = 'axis', text_args={'family':'sans-serif', 'color': 'white', 'ha': 'center', 'va': 'center', 'size':60, 'alpha':0.3, 'fontweight':'ultralight'})

    p.set_zlim(field = ('stars','particle_mass'), zmin = 5.e35, zmax = 1.e42)
    plot = p.plots[('stars','particle_mass')]
    plot.figure = fig
    plot.axes = grid[1].axes
    p._setup_plots()




    fig.set_size_inches(12, 6)
    fig.savefig('%s/cen_%s_%.4i.png'%(figdir, simname, DD))











































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
















