import numpy as np
from numpy import *
import matplotlib
import matplotlib.pyplot as plt

import yt
plt.ioff()
plt.close('all')

simnames = ['natural',
            'natural_v2',
            'natural_v3',
            'natural_v4',
            'nref11n_nref10f',
            'nref11c_nref9f',
            ]
clrs = ['black', 'grey', 'grey', 'grey', 'blue', 'red']
for DD in arange(600, 900, 50):
    DDname = 'DD%.4i'%DD
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for s, simname in enumerate(simnames):
        print (simname, DD)
        for aa, axs in enumerate(['x', 'y', 'z']):
            for i in np.arange(2):
                plunge = np.load('/Users/rsimons/Dropbox/rcs_foggie/freefall_experiments/plunge/%s_%s_%i_%s.npy'%(DDname, axs,i, simname), allow_pickle = True)[()]
                vmax = np.load('/Users/rsimons/Dropbox/rcs_foggie/freefall_experiments/vescape/%s_%s_vescape.npy'%(DDname, simname), allow_pickle = True)[()]
                dinner = yt.YTArray(200., 'kpc')
                dt = yt.YTArray(2.e7, 'yr')
                M = 0
                tot_Ms = []
                ts = []
                for t in arange(0, 1000):
                    douter = dinner
                    vmax_interp = yt.YTArray(np.interp(douter, vmax['r'], vmax['v']), 'km/s')
                    dinner = douter - (vmax_interp * dt.to('s')).to('kpc')
                    if dinner <0 : break
                    gd = where((plunge['d'] > dinner) & (plunge['d'] < douter))[0]
                    dvel = np.mean(plunge['vel'] + vmax_interp)
                    dens = np.mean(plunge['dens'])
                    P = dens * dvel**2.
                    M += P * dt
                    tot_Ms.append(M.value)
                    ts.append((t * dt.to('s')).to('Gyr'))
                ax.plot(ts, tot_Ms ,'-', color = clrs[s])    
                ax.set_yscale('log')
                ax.set_ylim(1.e-15, 1.e-9)
                ax.set_ylabel('Total Work')
                ax.set_xlabel('Time (Gyr)')

        #ax.annotate(simname, (0.1, 0.9), xycoords = 'axes fraction')
    ax.annotate(DDname, (0.05, 0.9), fontsize = 20, xycoords = 'axes fraction')

    fig.savefig('/Users/rsimons/Dropbox/rcs_foggie/freefall_experiments/figures/%s.png'%(DDname), dpi = 300)

