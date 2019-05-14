import os
import numpy as np
from numpy import *


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

for sat in sats:
    for a, ax in enumerate(axs):
        for z, zoom in enumerate(zooms):
            for d, DD in enumerate(DDs):
                for s, simname in enumerate(simnames):
                    fname = '/nobackupp2/rcsimons/foggie_momentum/sat_figures/%s/%s/%s/%s_%.4i_%.2i_%s_%s.png'%(simname, ax, zoom, simname, DD, sat, ax, zoom)
                    print os.path.isfile(fname)