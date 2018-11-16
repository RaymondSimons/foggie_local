import astropy
from astropy.io import fits
import numpy
from numpy import *
import matplotlib as mpl
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import yt
from scipy.signal import find_peaks
mpl.rcParams['text.usetex'] = True
mpl.rcParams['axes.linewidth'] = 5
plt.ioff()



def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = numpy.average(values, weights=weights)
    # Fast and numerically precise:
    variance = numpy.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))


def make_jzjcirc_fig(rd_s, ep_s, rd_d, ep_d, xmn, xmx, jmn, jmx, ms_s, ms_d, figname,   in_sel = None, alp1 = 1.0, alp2 = 1.0, bins_2d_stars = 1000, bins_2d_darkmatter = 50):
    fig, axes = plt.subplots(2,2, figsize = (25,15))


    axes[0,1].hist2d(rd_s, ep_s, range = ([xmn, xmx], [jmn, jmx]), norm = mpl.colors.LogNorm(), weights = ms_s, cmap = mpl.cm.viridis, bins = bins_2d_stars,      alpha = alp1)
    axes[1,1].hist2d(rd_d, ep_d, range = ([xmn, xmx], [jmn, jmx]), norm = mpl.colors.LogNorm(), weights = ms_d, cmap = mpl.cm.Greys,   bins = bins_2d_darkmatter, alpha = alp2)

    if in_sel is not None:
        axes[0,1].hist2d(rd_s[in_sel], ep_s[in_sel], range = ([xmn, xmx], [jmn, jmx]), norm = mpl.colors.LogNorm(), weights = ms_s[in_sel], cmap = mpl.cm.viridis, bins = bins_2d_stars, alpha = alp2)





    xmx = 25
    axes[0,0].hist2d(rd_s, ep_s, range = ([xmn, xmx], [jmn, jmx]), norm = mpl.colors.LogNorm(), weights = ms_s, cmap = mpl.cm.viridis, bins = bins_2d_stars,      alpha = alp1)
    axes[1,0].hist2d(rd_d, ep_d, range = ([xmn, xmx], [jmn, jmx]), norm = mpl.colors.LogNorm(), weights = ms_d, cmap = mpl.cm.Greys,   bins = bins_2d_darkmatter, alpha = alp2)

    if in_sel is not None:
        if rd_s[in_sel].min() < 25:
            axes[0,0].hist2d(rd_s[in_sel], ep_s[in_sel], range = ([xmn, xmx], [jmn, jmx]), norm = mpl.colors.LogNorm(), weights = ms_s[in_sel], cmap = mpl.cm.viridis, bins = bins_2d_stars, alpha = alp2)



    axes[0,0].annotate('stars', (0.80, 0.05), xycoords = 'axes fraction', fontsize = 45, fontweight = 'bold')
    axes[1,0].annotate('dark matter', (0.60, 0.05), xycoords = 'axes fraction', fontsize = 45, fontweight = 'bold')




    for ax in axes.ravel(): 
        ax.tick_params(axis='both', which='major', labelsize=20, size = 7, width = 3)
        ax.set_yticks([-2, -1, 0, 1, 2])


    for ax in axes[:,0]:    ax.set_ylabel(r'$\frac{j_{z}}{j_{circ}}$', fontsize = 50, labelpad = 40, rotation = 0)
    for ax in axes[1,:]:    ax.set_xlabel('distance from center (kpc)', fontsize = 30)
    for ax in axes[:,1]:    ax.set_ylabel('')
    for ax in axes[0,:]:    ax.set_xlabel('')


    axes[0,0].annotate(simname.replace('_', '-') +'\n'+'DD%.4i'%DDnum, (0.4, 0.9), xycoords = 'figure fraction',  fontsize = 40)

    fig.subplots_adjust(wspace = 0.10, hspace = 0.10)
    fig.savefig(out_fig_direct + '/' + figname, dpi = 400)


def make_2dhist_fig(rd_s, xmn, xmx, ms_s, figname, dist_mx_kpc = 1., thresh_mass_msun = 1.e5, thresh_rad_kpc = 0.2, ignore_rad_kpc = 5):
    fig, ax = plt.subplots(1,1, figsize = (18,9))

    hst, edges = histogram(rd_s, weights = ms_s, bins = linspace(0, xmx, bns))

    rds = np.mean([edges[0:len(edges) - 1], edges[1:len(edges)]], axis = 0)


    ax.plot(rds, hst, zorder = 2)


    peaks, _ = find_peaks(hst, distance = dist_mx_kpc*bns/(xmx-xmn), threshold = thresh_mass_msun/(thresh_rad_kpc*bns/(xmx-xmn)))


    peaks = peaks[rds[peaks] > ignore_rad_kpc]


    for p in peaks:
        ax.axvline(rds[p], color = 'black', zorder = 1, linestyle = '--', linewidth = 1, alpha = 0.6)
        
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel('distance from center (kpc)', fontsize = 30)


    ax.annotate(simname.replace('_', '-') +'\n'+'DD%.4i'%DDnum, (0.15, 0.9), xycoords = 'figure fraction',  fontsize = 30)


    ax.annotate('%i peaks detected'%len(peaks), (0.7, 0.93), xycoords = 'figure fraction',  fontsize = 25)


    ax.set_xlim(0, xmx)


    ax.axvspan(xmin = 0, xmax = ignore_rad_kpc, ymin = 0, ymax = 1, alpha = 0.2, color = 'grey')
    ax.annotate('peaks ignored', (ignore_rad_kpc*0.4 /(ax.get_xlim()[1] - ax.get_xlim()[0]), 0.3), xycoords = 'axes fraction',rotation = 90., color = 'black', fontweight = 'bold', fontsize = 20)


    fig.savefig(out_fig_direct + '/' + figname, dpi = 400)

    return peaks, rds, hst
        
        




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




    #make_jzjcirc_fig(rd_s, ep_s, rd_d, ep_d, xmn, xmx, jmn, jmx, ms_s, ms_d, figname = 'jzjcirc_stars_DM_DD%.4i.png'%DDnum, alp = 1.0)



    peaks, rds, hst = make_2dhist_fig(rd_s, xmn, xmx, ms_s, figname = 'stars_hist_rad_DD%.4i.png'%DDnum)




    for i,p in enumerate(peaks):
        mn_r = rds[p] - 0.5
        mx_r = rds[p] + 0.5
        in_sel = (rd_s > mn_r) & (rd_s < mx_r)
        #make_jzjcirc_fig(rd_s, ep_s, rd_d, ep_d, xmn, xmx, jmn, jmx, ms_s, ms_d, figname = 'sat_%.2i_jzjcirc_DD%.4i.png'%(i, DDnum),  in_sel = in_sel, alp1 = 0.2, alp2 = 1.0)


        x_s_avg, x_s_std = weighted_avg_and_std(x_s[in_sel], weights = ms_s[in_sel])
        y_s_avg, y_s_std = weighted_avg_and_std(y_s[in_sel], weights = ms_s[in_sel])
        z_s_avg, z_s_std = weighted_avg_and_std(z_s[in_sel], weights = ms_s[in_sel])


        print x_s_avg, x_s_std
        print y_s_avg, y_s_std
        print z_s_avg, z_s_std



        good       = where((abs(x_s[in_sel] - x_s_avg) < 3 * x_s_std) & 
                           (abs(y_s[in_sel] - y_s_avg) < 3 * y_s_std) & 
                           (abs(z_s[in_sel] - z_s_avg) < 3 * z_s_std))[0]


        good_ids  = id_s[in_sel][good]
        good_ages = age_s[in_sel][good]
        good_rds = rd_s[in_sel][good]
        good_eps = ep_s[in_sel][good]
        good_mss = ms_s[in_sel][good]
        good_xs = x_s[in_sel][good]
        good_ys = y_s[in_sel][good]
        good_zs = z_s[in_sel][good]
        good_vxs = vx_s[in_sel][good]
        good_vys = vy_s[in_sel][good]
        good_vzs = vz_s[in_sel][good]

        good_xs_box = x_s_box[in_sel][good]
        good_ys_box = y_s_box[in_sel][good]
        good_zs_box = z_s_box[in_sel][good]
        good_vxs_box = vx_s_box[in_sel][good]
        good_vys_box = vy_s_box[in_sel][good]
        good_vzs_box = vz_s_box[in_sel][good]








    '''
    for i,p in enumerate(peaks):



        axes[i].hist2d(rd_s, ep_s, range = ([xmn, xmx], [jmn,jmx]), norm = mpl.colors.LogNorm(), cmap = mpl.cm.viridis, bins = bns_2d, alpha = 0.2)

        axes[i].hist2d(rd_s[in_sel], ep_s[in_sel], range = ([xmn, xmx], [jmn,jmx]), norm = mpl.colors.LogNorm(), cmap = mpl.cm.viridis, bins = bns_2d)



        axes[i].annotate('stars', (0.8, 0.9), xycoords = 'axes fraction', fontsize = 50, fontweight = 'bold')
        axes[i].tick_params(axis='both', which='major', labelsize=20)
        axes[i].set_xlabel('distance from center (kpc)', fontsize = 30)
        axes[i].set_ylabel(r'$\frac{j_{z}}{j_{circ}}$', fontsize = 50, labelpad = 40, rotation = 0)

    '''






