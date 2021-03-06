import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
import astropy
from astropy.io import fits as pyfits
import glob
from glob import glob
import astrodendro
from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel, convolve_fft
import photutils 
from photutils import detect_sources
from photutils import *
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from joblib import Parallel, delayed
from astropy.io import fits
from numpy import *
import matplotlib as mpl
import os, sys, argparse


plt.ioff()
#This file will be used to store the profile of the momentum

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

    parser.add_argument('-run_parallel', '--run_parallel', default=False, help='Run parallel')

    parser.add_argument('-simname', '--simname', default=None, help='Simulation to be analyzed.')

    parser.add_argument('-snapname', '--snapname', default=None, help='Snapshot files to be analyzed.')

    parser.add_argument('-haloname', '--haloname', default='halo_008508', help='halo_name')
    parser.add_argument('-on_system', '--on_system', default='pfe', help='System being used (pfe or local)')
    parser.add_argument('-ddmin', '--ddmin', default=500, help='halo_name')
    parser.add_argument('-ddmax', '--ddmax', default=600, help='halo_name')
    parser.add_argument('-n_jobs', '--n_jobs', default=3, help='number of jobs')



    args = vars(parser.parse_args())
    return args



def write_fits(fits_name, mom_data, merger_tag, x_stars_box , y_stars_box , z_stars_box, vx_stars_box , vy_stars_box , vz_stars_box):

    print '\tGenerating fits for %s'%fits_name
    master_hdulist = []
    master_hdulist.append(mom_data['PRIMARY'])
    colhdr = fits.Header()
    master_hdulist.append(mom_data['nir_mstar_cat'])
    master_hdulist.append(mom_data['nir_net_momentum'])
    master_hdulist.append(mom_data['nir_net_momentum_s'])
    master_hdulist.append(mom_data['stars_id'])
    master_hdulist.append(fits.ImageHDU(data = np.stack((x_stars_box , y_stars_box , z_stars_box,)), header = colhdr, name = 'stars_xyz_box_position'))
    master_hdulist.append(fits.ImageHDU(data = np.stack((vx_stars_box , vy_stars_box , vz_stars_box)), header = colhdr, name = 'stars_xyz_box_velocity'))
    master_hdulist.append(mom_data['star_mass'])
    master_hdulist.append(mom_data['star_age'])
    master_hdulist.append(fits.ImageHDU(data = merger_tag, header = colhdr, name = 'star_merger_tag'))

    print '\tSaving to ' + fits_name
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto(fits_name, clobber = True)

    return master_hdulist


def make_heatmap(ax, epsilon, zz_gas, min_z, max_z, weights = None, good = None, xlabel = 'z height (kpc)', ylabel = 'j$_z$/j$_{circ}$', bins_n = 200, eps_min = 2, eps_max = 2, segm = None, srt_labels = [], do_plot = True):
    #if weights == None:
    #    weights = np.ones(len(zz_gas))

    if good:        
        epsilon = epsilon[good]
        zz_gas = zz_gas[good]
        weights = weights[good]
    heatmap, xedges, yedges = np.histogram2d(epsilon, zz_gas, bins=[linspace(eps_min,eps_max,bins_n), linspace(min_z,max_z,bins_n)], weights = weights)




    srt_labels = array(srt_labels)

    heatmap_sorted = sort(heatmap.ravel())
    heatmap_sorted = heatmap_sorted[heatmap_sorted > 0.]

    vmn_scale = 0.1
    vmx_scale = 0.998
    vmn = heatmap_sorted[int(vmn_scale*len(heatmap_sorted))]
    vmx = heatmap_sorted[int(vmx_scale*len(heatmap_sorted))]
    #heatmap = np.ma.masked_where((heatmap < 1), heatmap)
    #heatmap.data[heatmap.data < 1.] = nan

    #heatmap.data[segm > 1] = 0
    if len(srt_labels)>0:
        #for lbl in srt_labels[1:len(srt_labels)]:
        #    heatmap.data[segm == lbl] = 0
        heatmap.data[segm!=srt_labels[0]] = 0

    if do_plot:
        ax.imshow(heatmap, interpolation = 'nearest', norm = mpl.colors.LogNorm(vmin = vmn, vmax = vmx), origin = 'lower', cmap = 'viridis')
        kern = Gaussian2DKernel(1.)
        kern.normalize()
        heatmap_conv = convolve_fft(heatmap, kern)
        heatmap_conv = np.ma.masked_where((heatmap_conv < 10), heatmap_conv)
        heatmap_conv.data[heatmap_conv.data < 10.] = nan

        if False:
            X = arange(heatmap.data.shape[0])
            Y = arange(heatmap.data.shape[1])
            Z = log10(heatmap.data)
            ax.contour(X, Y, Z, 3, colors = 'grey')



        ax.set_yticks([0,(bins_n-1)/4,(bins_n-1)/2,3*(bins_n-1)/4,bins_n-1])
        #ax.set_xticks([0,(bins_n-1)/2,bins_n-1])
        #ax.set_xticklabels([format(yedges[0],'.0f'),format(yedges[(bins_n-1)/2],'.0f'),format(yedges[bins_n-1],'.0f')])

        xtcks = [min_z, max_z/2., max_z]
        xtcks_c = [xt * bins_n/(1.*max_z - min_z) for xt in xtcks]
        xtcks_l = [format(xt, '.0f') for xt in xtcks]


        ytcks = np.arange(eps_min, eps_max + 1, 1)
        ytcks_c = [(yt - eps_min) * bins_n/(1.*eps_max - eps_min) for yt in ytcks]
        ytcks_l = [format(yt, '.0f') for yt in ytcks]

        ax.set_xticks(xtcks_c)
        ax.set_xticklabels(xtcks_l)

        print ytcks, ytcks_c, ytcks_l
        ax.set_yticks(ytcks_c)
        ax.set_yticklabels(ytcks_l)


        #ax.set_yticklabels([''])
        #ax.set_yticklabels([format(xedges[0],'.0f'),format(xedges[int((bins_n-1)/4)],'.0f'), format(xedges[int((bins_n-1)/2.)],'.0f'),format(xedges[int(3*(bins_n-1)/4.)],'.0f'),format(xedges[bins_n-1],'.0f')])
        #ax.set_xticklabels([''])
        ax.set_xlabel(xlabel, fontsize = 15)

        ax.set_ylabel(ylabel, fontsize = 20)
        ax.minorticks_on()

        ax.tick_params(axis="both", which='major', color='black', labelcolor='black',size=5, width=1.5)
        ax.tick_params(axis="both", which='minor', color='black', labelcolor='black',size=3, width=1.5)

        return ax, heatmap
    else:
        return heatmap



def add_at(ax, t, loc=2):
    fp = dict(size=10)
    _at = AnchoredText(t, loc=loc, prop=fp)
    ax.add_artist(_at)
    return _at


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    values = values[~isnan(weights)]
    weights = weights[~isnan(weights)]
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))


def find_thresh(mn, mx, npix, heatmap):
    nlabels = 0.
    segm_labels_prev = 0
    mr_prev2 = -99
    mr_prev = -99
    kern = Gaussian2DKernel(0.2, x_size = 4*10, y_size = 4*10)
    kern.normalize()
    a = zeros(kern.array.shape)
    a[int(kern.array.shape[1]/2.),int(kern.array.shape[1]/2.)] = 1
    kern_2 = Gaussian1DKernel(8)
    a[:,int(kern.array.shape[1]/2.)] = convolve_fft(a[:,int(kern.array.shape[1]/2.)], kern_2)
    a/=sum(a)
    b = convolve_fft(a, kern)
    b/=sum(b)
    temp_heatmap = convolve_fft(heatmap.data, b)
    temp_heatmap[temp_heatmap <= 0] = nan

    for tt, t in enumerate(linspace(mn, mx, 1000)):

        threshold = t
        segm = detect_sources(log10(temp_heatmap), threshold = threshold, npixels = npix)  
        masses = array([sum(temp_heatmap[segm.array == lbl]) for lbl in arange(1, segm.nlabels+1)])
        srt_masses = masses[argsort(masses)[::-1]]
        if len(masses) > 1:
            mass_ratio = srt_masses[0]/srt_masses[1]
            if mr_prev == -99:
                mr_prev = mass_ratio  
                thresh = threshold       
            if (log10(srt_masses[0]) > 7.5) & (log10(srt_masses[1]) > 7.5) & \
                (mr_prev/mass_ratio > 10) & (mass_ratio < 100) & (nansum(srt_masses) > 0.50*nansum(temp_heatmap)):
                thresh = threshold


            mr_prev = mass_ratio

            if len(masses) > 2:
                mass_ratio2 = srt_masses[0]/srt_masses[2]
                if mr_prev2 == -99:
                    mr_prev2 = mass_ratio2
                    thresh = threshold    
                    
                if (log10(srt_masses[0]) > 7.5) & (log10(srt_masses[1]) > 7.5) & (mr_prev2/mass_ratio2 > 10) & (mass_ratio2 < 300) & (nansum(srt_masses) > 0.50*nansum(temp_heatmap)):
                    thresh = threshold
                    
                mr_prev2 = mass_ratio2
        segm_labels_prev = segm.nlabels
    return thresh, temp_heatmap


def make_figure(snap_name, simname, on_system = 'pfe'):
    if True:
        if on_system == 'pfe':
            mom_outdir = '/nobackupp2/rcsimons/foggie_momentum/momentum_fits'
            fig_file = '/nobackupp2/rcsimons/foggie_momentum/figures/%s_%s_satellites.png'%(simname, snap_name)
        else:
            mom_outdir = '/Users/rsimons/Dropbox/rcs_foggie/outputs'
            fig_file = '/Users/rsimons/Dropbox/rcs_foggie/outputs/%s_%s_satellites.png'%(simname, snap_name)

        mom_fl = '%s/%s_%s_momentum.fits'%(mom_outdir, simname, snap_name)
        mom_data = fits.open(mom_fl)

        eps_min = -2.
        eps_max = 2.
        r_min = 0.
        r_max = 30
        bins_n_stars = 2000
        bins_n_dark = 30

        max_nmergers = 20


        x_pos_stars, y_pos_stars, z_pos_stars = mom_data['STARS_GAL_POSITION'].data
        x_vel_stars, y_vel_stars, z_vel_stars = mom_data['STARS_GAL_VELOCITY'].data
        r_pos_stars = sqrt(x_pos_stars**2. + y_pos_stars**2. + z_pos_stars**2.)
        epsilon_stars = mom_data['STARS_EPSILON_FIXED'].data


        x_pos_dark, y_pos_dark, z_pos_dark = mom_data['DARK_GAL_POSITION'].data
        x_vel_dark, y_vel_dark, z_vel_dark = mom_data['DARK_GAL_VELOCITY'].data
        r_pos_dark = sqrt(x_pos_dark**2. + y_pos_dark**2. + z_pos_dark**2.)
        epsilon_dark = mom_data['DARK_EPSILON_FIXED'].data



        '''
        epsilon_stars_digitized = np.digitize(epsilon_stars, bins = linspace(eps_min, eps_max, bins_n))
        r_stars_digitized = np.digitize(r_pos_stars, bins = linspace(r_min, r_max, bins_n))
        empt_arr = np.empty((bins_n-1,bins_n-1), dtype = object)

        for i in arange(bins_n-1):
            good_r_stars = where(r_stars_digitized == i)[0]
            r_stars_digitized_new = r_stars_digitized[good_r_stars]
            epsilon_stars_digitized_new = epsilon_stars_digitized[good_r_stars]                
            for j in arange(bins_n-1):
                good_eps_stars = good_r_stars[where(epsilon_stars_digitized_new == j)[0]]
                empt_arr[i,j] = good_eps_stars
        '''

        star_age = mom_data['STAR_AGE'].data
        star_mass= mom_data['STAR_MASS'].data

        dark_age = mom_data['DARK_AGE'].data
        dark_mass= mom_data['DARK_MASS'].data




    if True:
        fig, axes  = plt.subplots(2,2, figsize = (10, 10))


        good = where(star_age < 2.e40)





        axes[0,0], heatmap   = make_heatmap(axes[0,0], epsilon_stars, r_pos_stars, min_z = r_min, max_z = r_max, weights = star_mass, 
                            good = good, xlabel = '', ylabel = '', bins_n = bins_n_stars, eps_min = eps_min, eps_max = eps_max)


        axes[0,1], heatmap   = make_heatmap(axes[0,1], epsilon_stars, r_pos_stars, min_z = r_min, max_z = 100., weights = star_mass, 
                            good = good, xlabel = '', ylabel = '', bins_n = bins_n_stars, eps_min = eps_min, eps_max = eps_max)

        good = where(dark_age < 2.e40)

        axes[1,0], heatmap   = make_heatmap(axes[1,0], epsilon_dark, r_pos_dark, min_z = r_min, max_z = r_max, weights = dark_mass, 
                            good = good, xlabel = '', ylabel = '', bins_n = bins_n_dark, eps_min = eps_min, eps_max = eps_max)


        axes[1,1], heatmap   = make_heatmap(axes[1,1], epsilon_dark, r_pos_dark, min_z = r_min, max_z = 100., weights = dark_mass, 
                            good = good, xlabel = '', ylabel = '', bins_n = bins_n_dark, eps_min = eps_min, eps_max = eps_max)

    for ax in axes.ravel():
        ax.set_xlabel('radial distance (kpc)', fontsize = 18)
        ax.set_ylabel(r'$j_z/j_{circ}$', fontsize = 18)




    fig.tight_layout() 

    plt.savefig(fig_file, dpi = 300)



if __name__ == '__main__':    
        #simname = 'nref11n_selfshield_z15'

        args = parse()
        simname =  args['simname']
        snap_name =  args['snapname']
        on_system = args['on_system']
        #Parallel(n_jobs = 5, backend = 'threading')(delayed(make_figure)(snap_name = 'DD%.4i'%i, simname = simname) for i in np.arange(500, 505))
        make_figure(snap_name = snap_name, simname = simname, on_system = on_system)























