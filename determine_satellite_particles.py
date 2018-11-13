import astropy
from astropy.io import fits
import numpy
from numpy import *
import matplotlib as mpl
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import yt
from scipy.signal import find_peaks
plt.ioff()




def make_jzjcirc_fig(rd_s, ep_s, rd_d, ep_d, xmn, xmx, jmn, jmx, ms_s, ms_d, figname,  alp = 1.0, bins_2d_stars = 1000, bins_2d_darkmatter = 50):
    fig, axes = plt.subplots(2,2, figsize = (25,15))


    axes[0,1].hist2d(rd_s, ep_s, range = ([xmn, xmx], [jmn, jmx]), norm = mpl.colors.LogNorm(), weights = ms_s, cmap = mpl.cm.viridis, bins = bins_2d_stars,      alpha = alp)
    axes[1,1].hist2d(rd_d, ep_d, range = ([xmn, xmx], [jmn, jmx]), norm = mpl.colors.LogNorm(), weights = ms_d, cmap = mpl.cm.Greys,   bins = bins_2d_darkmatter, alpha = alp)


    xmx = 25
    axes[0,0].hist2d(rd_s, ep_s, range = ([xmn, xmx], [jmn, jmx]), norm = mpl.colors.LogNorm(), weights = ms_s, cmap = mpl.cm.viridis, bins = bins_2d_stars,      alpha = alp)
    axes[1,0].hist2d(rd_d, ep_d, range = ([xmn, xmx], [jmn, jmx]), norm = mpl.colors.LogNorm(), weights = ms_d, cmap = mpl.cm.Greys,   bins = bins_2d_darkmatter, alpha = alp)




    axes[0,0].annotate('stars', (0.75, 0.05), xycoords = 'axes fraction', fontsize = 40, fontweight = 'bold')
    axes[1,0].annotate('dark matter', (0.55, 0.05), xycoords = 'axes fraction', fontsize = 40, fontweight = 'bold')




    for ax in axes.ravel(): 
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.set_yticks([-2, -1, 0, 1, 2])


    for ax in axes[:,0]:    ax.set_ylabel(r'$\frac{j_{z}}{j_{circ}}$', fontsize = 50, labelpad = 40, rotation = 0)
    for ax in axes[1,:]:    ax.set_xlabel('distance from center (kpc)', fontsize = 30)
    for ax in axes[:,1]:    ax.set_ylabel('')
    for ax in axes[0,:]:    ax.set_xlabel('')


    axes[0,0].annotate('%s\nDD%.4i'%(simname, DDnum), (0.4, 0.9), xycoords = 'figure fraction',  fontsize = 40)

    fig.subplots_adjust(wspace = 0.10, hspace = 0.10)
    fig.savefig(out_fig_direct + '/' + figname, dpi = 400)






if __name__ == '__main__':

    xmn, xmx = 0, 150
    jmn, jmx = -2, 2


    global out_fig_direct, DDnum, simname
    simname = 'nref11n_selfshield_z15'
    out_fig_direct = '/Users/rsimons/Dropbox/rcs_foggie/sat_figures/%s'%simname
    DDnum = 350

    bns = 500

    a = fits.open('/Users/rsimons/Dropbox/rcs_foggie/outputs/momentum_fits/%s_DD%.4i_momentum.fits'%(simname, DDnum))


    ms_s = a['STAR_MASS'].data
    age_s = a['STAR_AGE'].data

    id_s = a['STARS_ID'].data
    rd_s = sqrt(sum(a['STARS_GAL_POSITION'].data**2., axis = 0))
    x_s, y_s, z_s = a['STARS_GAL_POSITION'].data
    vx_s, vy_s, vz_s = a['STARS_GAL_VELOCITY'].data

    x_s_box, y_s_box, z_s_box = a['STARS_BOX_POSITION'].data
    vx_s_box, vy_s_box, vz_s_box = a['STARS_BOX_VELOCITY'].data



    ep_s = a['STARS_EPSILON_FIXED'].data

    ms_d = a['DARK_MASS'].data
    rd_d = sqrt(sum(a['DARK_GAL_POSITION'].data**2., axis = 0))
    x_d, y_d, z_d = a['DARK_GAL_POSITION'].data
    vx_d, vy_d, vz_d = a['DARK_GAL_VELOCITY'].data

    ep_d = a['DARK_EPSILON_FIXED'].data




    make_jzjcirc_fig(rd_s, ep_s, rd_d, ep_d, xmn, xmx, jmn, jmx, ms_s, ms_d, figname = 'jzjcirc_stars_DM_DD%.4i.png'%DDnum, alp = 1.0)











