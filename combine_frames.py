import os
import numpy as np
from numpy import *
from PIL import Image


DDs = arange(49, 1000)
simnames =  ['natural',
            'natural_v2',
            'natural_v3',
            'natural_v4',
            'nref11c_nref9f',
            'nref11n_nref10f']
axs = ['x', 'y', 'z']
zooms = ['zoomin', 'zoomout', 'zoomoutfar']
sats = arange(6)



DDs = arange(49, 50)
simnames =  ['natural',
            'natural_v2',
            'natural_v3',
            'natural_v4',
            'nref11c_nref9f',
            'nref11n_nref10f']
axs = ['x']
zooms = ['zoomin']
sats = arange(1)


for sat in sats:
    for a, ax in enumerate(axs):
        for z, zoom in enumerate(zooms):
            for d, DD in enumerate(DDs):
                imgs = [] 
                fname_out = '/nobackupp2/rcsimons/foggie_momentum/sat_figures/combined/%s/%s/%.4i_%.2i_%s_%s.png'%(ax, zoom, DD, sat, ax, zoom)
                for s, simname in enumerate(simnames):
                    fname = '/nobackupp2/rcsimons/foggie_momentum/sat_figures/%s/%s/%s/%s_%.4i_%.2i_%s_%s.png'%(simname, ax, zoom, simname, DD, sat, ax, zoom)
                    print (os.path.isfile(fname))
                    im =  Image.open(fname)
                    imgs.append(im)
                #imgs_12 = imgs[0:2]
                #imgs_34 = imgs[2:4]
                #imgs_56 = imgs[4:6]
                #imgs_comb = []

                #for imgs in [imgs_12, imgs_34, imgs_56]:
                min_shape = sorted([(np.sum(i.size), i.size ) for i in imgs])[0][1]
                imgs_comb_temp = np.vstack((np.asarray( i.resize(min_shape) ) for i in imgs ))
                imgs_comb_temp = Image.fromarray( imgs_comb)
                #imgs_comb.append(imgs_comb_temp)
                imgs_comb_temp.save(fname_out)

                '''
                min_shape = sorted([(np.sum(i.size), i.size ) for i in imgs_comb])[0][1]
                imgs_comb = np.hstack((np.asarray( i.resize(min_shape) ) for i in imgs_comb ))
                imgs_comb_all = Image.fromarray( imgs_comb)
                imgs_comb_all.save(fname_out)
    
                '''



















