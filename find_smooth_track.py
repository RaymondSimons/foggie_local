import matplotlib
import astropy
from astropy.io import fits
import numpy as np
from numpy import *
import math
from joblib import Parallel, delayed
import os, sys, argparse
import yt
import matplotlib.pyplot as plt

plt.ioff()
plt.close('all')


gp_dir = '/Users/rsimons/Dropbox/rcs_foggie/galprops/halo_008508'

simnames = ['nref11n_selfshield_z15', 'nref11n_nref10f_selfshield_z6']



for simname in simnames:
    fig, axes = plt.subplots(1,3, figsize = (12, 4))
    x = []
    y = []
    z = []
    DD_arr = arange(200, 1056)
    for DD in DD_arr:
        try:
            a = np.load(gp_dir + '/%s_DD%.4i_galprops.npy'%(simname, DD))[()]
            x.append(a['stars_center'][0][0])
            y.append(a['stars_center'][0][1])
            z.append(a['stars_center'][0][2])
        except:

            x.append(nan)
            y.append(nan)
            z.append(nan)

    x = array(x)
    y = array(y)
    z = array(z)

    DD_arr = DD_arr[~isnan(z)]
    x = x[~isnan(z)]
    y = y[~isnan(z)]
    z = z[~isnan(z)]

    x_fit = polyfit(DD_arr, x, deg = 4)
    y_fit = polyfit(DD_arr, y, deg = 4)
    z_fit = polyfit(DD_arr, z, deg = 4)


    axes[0].plot(DD_arr, x, 'k', linewidth = 1)
    axes[1].plot(DD_arr, y, 'k', linewidth = 1)
    axes[2].plot(DD_arr, z, 'k', linewidth = 1)


    xf = x_fit[0]*DD_arr**4. +x_fit[1]*DD_arr**3. + x_fit[2]*DD_arr**2. + x_fit[3]*DD_arr**1. + x_fit[4]
    yf = y_fit[0]*DD_arr**4. +y_fit[1]*DD_arr**3. + y_fit[2]*DD_arr**2. + y_fit[3]*DD_arr**1. + y_fit[4]
    zf = z_fit[0]*DD_arr**4. +z_fit[1]*DD_arr**3. + z_fit[2]*DD_arr**2. + z_fit[3]*DD_arr**1. + z_fit[4]



    axes[0].plot(DD_arr, xf, 'b--', linewidth = 0.4)
    axes[1].plot(DD_arr, yf, 'b--', linewidth = 0.4)
    axes[2].plot(DD_arr, zf, 'b--', linewidth = 0.4)


    fit_params = {}
    fit_params['x'] = x_fit
    fit_params['y'] = y_fit
    fit_params['z'] = z_fit
    fit_params['DDarr'] = DD_arr

    np.save('/Users/rsimons/Dropbox/rcs_foggie/catalogs/center_%s.npy'%simname, fit_params)



    fig.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/%s_trackcen.png'%simname, dpi = 400)













    