import astropy
from astropy.io import fits
import yt
from yt import YTArray
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from numpy import *
from consistency import *

plt.ioff()
plt.close('all')

ds = yt.load('/user/rsimons/foggie/sims/halo_008508/nref11n_nref10f/DD0809/DD0809')



wd = 30.
a = np.load('/user/rsimons/foggie_outputs/nref11n_nref10f_interpolations_DD0150_new.npy', allow_pickle = True)[()]

sat = 1
cen_x = a['SAT_%.2i'%sat]['fx'](809)  
cen_y = a['SAT_%.2i'%sat]['fy'](809)  
cen_z = a['SAT_%.2i'%sat]['fz'](809)  


cen_g = YTArray([cen_x, cen_y, cen_z], 'kpc')

x_w = wd
y_w = wd
z_w = wd

W = YTArray([x_w, y_w, z_w], 'kpc')



for i in arange(300, 301):
    print (i)
    box = ds.r[cen_g[0] - yt.YTArray(wd,  'kpc'): cen_g[0] + yt.YTArray( wd, 'kpc'), \
               cen_g[1] - yt.YTArray(wd, 'kpc'):  cen_g[1] +  yt.YTArray(wd, 'kpc'), \
               cen_g[2] - yt.YTArray(wd, 'kpc'):  cen_g[2] +  yt.YTArray(wd, 'kpc')]

    '''
    L = [1*cos(pi*(i)/100.),0, 1*sin(pi*(i)/100.)] # vector normal to cutting plane

    N = 512
    north_vector = [0,1,0]
    fig, ax = plt.subplots(1,1, figsize = (10, 10))
    ax.axis('off')


    image1 = yt.off_axis_projection(box, cen_g, L, W, N, ('gas', 'density'), north_vector =  north_vector)
    '''
    axis = 'x'
    p = yt.ProjectionPlot(ds, axis, ("gas","density"), center = cen_g, data_source=box, width=W)
    p.set_unit(('gas','density'), 'Msun/pc**2')
    p.set_zlim(('gas', 'density'), zmin = density_proj_min, zmax =  density_proj_max)
    p.set_cmap(('gas', 'density'), density_color_map)
    p.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
    p.hide_axes()
    p.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
    p.annotate_scale(size_bar_args={'color':'white'})
    p.save("/user/rsimons/foggie/figures/movies/ram_pressure/frame_%i.png"%i)
    '''
    image1 = image1.in_units('Msun/pc**2')
    ax.imshow(np.log10(image1), cmap = density_color_map)    
    fig.subplots_adjust(left = 0.0, right = 0.92, top =1.0, bottom = 0.0, hspace = 0.0, wspace = 0.0)
    fig.savefig("/user/rsimons/foggie/figures/movies/ram_pressure/frame_%i.png"%i, dpi = 300)
    '''













