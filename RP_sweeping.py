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



wd = 10.

x_w = wd
y_w = wd
z_w = wd

W = YTArray([x_w, y_w, z_w], 'kpc')

cen_g = YTArray([31597.918859155205, 31096.55501686099, 32192.674736309542], 'kpc')


for i in arange(300, 301):
    print (i)
    L = [1*cos(pi*(i)/100.),0, 1*sin(pi*(i)/100.)] # vector normal to cutting plane

    N = 512
    north_vector = [0,1,0]
    fig, ax = plt.subplots(1,1, figsize = (10, 10))
    ax.axis('off')

    box = ds.r[cen_g[0] - yt.YTArray(wd,  'kpc'): cen_g[0] + yt.YTArray( wd, 'kpc'), \
               cen_g[1] - yt.YTArray(wd, 'kpc'):  cen_g[1] +  yt.YTArray(wd, 'kpc'), \
               cen_g[2] - yt.YTArray(wd, 'kpc'):  cen_g[2] +  yt.YTArray(wd, 'kpc')]


    image1 = yt.off_axis_projection(box, cen_g, L, W, N, ('gas', 'density'), north_vector =  north_vector)


    image1 = image1.in_units('Msun/pc**2')
    ax.imshow(np.log10(image1), cmap = density_color_map)    
    fig.subplots_adjust(left = 0.0, right = 0.92, top =1.0, bottom = 0.0, hspace = 0.0, wspace = 0.0)
    fig.savefig("/user/rsimons/foggie/figures/movies/ram_pressure/frame_%i.png"%i, dpi = 300)














