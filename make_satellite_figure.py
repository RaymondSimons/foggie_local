import matplotlib
matplotlib.use('Agg')
import astropy
from astropy.io import fits
import numpy as np
from numpy import *
import math
from joblib import Parallel, delayed
import os, sys, argparse
import yt
import matplotlib.pyplot as plt
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
    parser.add_argument('-DD', '--DD', default=None, help='DD to use')

    parser.add_argument('-simname', '--simname', default=None, help='Simulation to be analyzed.')

    parser.add_argument('-haloname', '--haloname', default='halo_008508', help='halo_name')

    args = vars(parser.parse_args())
    return args



if __name__ == '__main__':
    args = parse()

    simname = args['simname']
    DD = int(args['DD'])
    haloname = args['haloname']

    snapname = 'DD%.4i'%DD

    ds = yt.load('/nobackupp2/mpeeples/%s/%s/%s/%s'%(haloname, simname, snapname, snapname))

    #cen_file =  np.load('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_cen.npy'%(simname, DD))[()]

    if simname == 'natural': cen_name = 'natural'
    if 'v2' in simname: cen_name = 'natural_v2'
    if 'v3' in simname: cen_name = 'natural_v3'
    if 'v4' in simname: cen_name = 'natural_v4'

    cen_fits = fits.open('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s/%s_DD%.4i_anchorprops.fits'%(cen_name, simname, DD))


    #anchor_xs_box_avg, anchor_ys_box_avg, anchor_zs_box_avg, anchor_vxs_box_avg, anchor_vys_box_avg, anchor_vzs_box_avg = cen_file
    cenx = cen_fits['SAT_%.2i'%sat_n].data['box_avg'][0]
    ceny = cen_fits['SAT_%.2i'%sat_n].data['box_avg'][1]
    cenz = cen_fits['SAT_%.2i'%sat_n].data['box_avg'][2]
    cen = yt.YTArray([cenx, ceny, cenz], 'kpc')

    #vel_vec = array([anchor_vxs_box_avg, anchor_vys_box_avg, anchor_vzs_box_avg])
    #L = vel_vec/np.linalg.norm(vel_vec)


    #L = [1*math.cos(pi*(DD)/100.),0, 1*math.sin(pi*(DD)/100.)] # vector normal to cutting plane

    Ls = [[1,0,0], [0,1,0], [0, 0, 1]]

    W1 = yt.YTArray([15, 15, 15], 'kpc')
    W2 = yt.YTArray([100, 100, 100], 'kpc')

    north_vector = [0,0.7,0.7]
    N = 512


    fig, axes = plt.subplots(3,2, figsize = (10.8, 15))



    for i in arange(3):
        L = Ls[i]
        image1 = yt.off_axis_projection(ds, cen, L, W1, N, ('gas', 'density'), north_vector =  north_vector)
        image2 = yt.off_axis_projection(ds, cen, L, W2, N, ('gas', 'density'), north_vector =  north_vector)

        image1 = image1.in_units('Msun * kpc**-2')
        image2 = image2.in_units('Msun * kpc**-2')

        im1 = axes[i, 0].imshow(np.log10(image1), vmin = 4.5, vmax = 9.5)
        im2 = axes[i, 1].imshow(np.log10(image2), vmin = 4.5, vmax = 9.5)


        bar_len_kpc = 1.
        bar_len_pix = 1.*N/W1[0].value * bar_len_kpc
        y_bar_start_pix = 0.5*N
        y_bar_end_pix = y_bar_start_pix + bar_len_pix
        x_bar_pix = 0.9*N
        axes[i, 0].plot([x_bar_pix, x_bar_pix], [y_bar_start_pix, y_bar_end_pix], color = 'white', linewidth = 3)
        axes[i, 0].annotate("%i kpc"%bar_len_kpc, (x_bar_pix-0.08*N, y_bar_start_pix- 0.04*N), color = 'white', fontsize = 15, rotation = 0)



        bar_len_kpc = 10
        bar_len_pix = 1.*N/W2[0].value * bar_len_kpc
        y_bar_start_pix = 0.5*N
        y_bar_end_pix = y_bar_start_pix + bar_len_pix
        x_bar_pix = 0.9*N
        axes[i, 1].plot([x_bar_pix, x_bar_pix], [y_bar_start_pix, y_bar_end_pix], color = 'white', linewidth = 3)
        axes[i, 1].annotate("%i kpc"%bar_len_kpc, (x_bar_pix-0.10*N, y_bar_start_pix- 0.04*N), color = 'white', fontsize = 15, rotation = 0)






    for ax in axes.ravel():
        ax.axis('off')

    cax = fig.add_axes([0.90, 0.0, 0.02, 1.0])
    cbr = fig.colorbar(im1, cax=cax,orientation="vertical")
    cbr.set_label('projected gas density (M$_{\odot}$ kpc$^{-2}$)', fontsize = 15)
    fig.subplots_adjust(left = 0.0, right = 0.92, top =1.0, bottom = 0.0, hspace = 0.0, wspace = 0.0)
    fig.savefig('/nobackupp2/rcsimons/foggie_momentum/sat_figures/%s_%s_0.png'%(simname, snapname), dpi = 500)



























