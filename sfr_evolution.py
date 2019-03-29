import astropy
from astropy.io import fits
import glob
from glob import glob
import scipy
import os
from scipy import interpolate
#plt.rcParams['text.usetex'] = True
plt.ioff()
plt.close('all')


selection_DD = 1049


#fits_dir = '/Users/rsimons/Dropbox/rcs_foggie/satellite_masses/selected_at_DD%s'%selection_DD
fits_dir = '/Users/rsimons/Dropbox/rcs_foggie/satellite_masses/all'
fig_dir = '/Users/rsimons/Dropbox/rcs_foggie/figures/sat_properties'
dist_from_center_file = '/Users/rsimons/Dropbox/rcs_foggie/catalogs/satellite_dist_from_center.npy'
dist_from_center = np.load(dist_from_center_file)[()]

DD_to_t = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/DD_time.npy')[()]
clrs = ['blue', 'darkblue', 'red']
prev = 0
#sat_ns_arr = [(0, 0), (4, 2), (3, 4)]
sat_ns_arr = [(4, 2)]

sat_ns_arr = [(x, x) for x in arange(11)]
 
sat_ns_arr = [(0, 0), (1, 1), 
              (2, 2), (3, 3), 
              (4, 4), (5, 5), 
              (6, 3), (7, 2), 
              (10,8), (11,7)]

sat_ns_arr = [(0, 0, 0), (1,1, 1), (2,2, 2), (3,3, 3), (4,4, 4), (5,5, 5), ]


for s, sat_ns in enumerate(sat_ns_arr):
    fig1, ax1 = plt.subplots(1,1, figsize = (6,6))
    fig2, ax2 =  plt.subplots(1,1, figsize = (6,6))
    fig3, ax3 =  plt.subplots(1,1, figsize = (6,6))
    fig4, ax4 =  plt.subplots(1,1, figsize = (6,6))


    #for sm, simname in enumerate(['natural', 'nref11n_nref10f_selfshield_z6']):
    for sm, simname in enumerate(['natural', 'nref11n_v2_selfshield_z15', 'nref11n_nref10f_selfshield_z6']):

        sat_n = sat_ns[sm]
        DD_arr_all = arange(44, 500)
        tt_arr_all = array([DD_to_t[2][where(DD_to_t[0] == Dnum)[0][0]] for Dnum in DD_arr_all])
        zz_arr_all = array([DD_to_t[1][where(DD_to_t[0] == Dnum)[0][0]] for Dnum in DD_arr_all])
        sfr_arr = []
        sms_arr = []
        gms_arr = []

        tt_arr = []
        zz_arr = []
        R_90 = []

        R_c = []
        R_cen = dist_from_center[max(sm-1, 0), :, sat_n]
        for d, DD in enumerate(DD_arr_all):
            #fits_fl = '%s_DD%.4i_mass_sat%.2i_%i.fits'%(simname, DD, sat_n, selection_DD)
            fits_fl = '%s_DD%.4i_mass_sat%.2i.fits'%(simname, DD, sat_n)
            if os.path.isfile(fits_dir + '/' + fits_fl):
                data = fits.open(fits_dir + '/' + fits_fl)
                #print data[6].data[-1]
                if data['STARS_MASS'].data[0] < 0.9 * data['STARS_MASS'].data[25]:
                    frac_mass = data['STARS_MASS'].data[:]/data['STARS_MASS'].data[25]
                    f = interpolate.interp1d(frac_mass, data['DISTANCE'].data)
                    f2 = interpolate.interp1d(data['DISTANCE'].data, data['STARS_YOUNGMASS'].data/(2.e7))
                    f3 = interpolate.interp1d(data['DISTANCE'].data, data['STARS_MASS'].data)
                    f4 = interpolate.interp1d(data['DISTANCE'].data, data['GAS_TOT'].data)

                    R_90.append(f(0.9))
                    #sfr_arr.append(f2(f(0.9)))
                    #sms_arr.append(f3(f(0.9)))
                    #gms_arr.append(f4(f(0.9)))


                    sfr_arr.append(data['STARS_YOUNGMASS'].data[30])
                    sms_arr.append(data['STARS_MASS'].data[30])
                    gms_arr.append(data['GAS_TOT'].data[30])



                    tt_arr.append(tt_arr_all[d])
                    zz_arr.append(zz_arr_all[d])

                    R_c.append(R_cen[d] * (1+zz_arr_all[d]))



        if simname == 'natural': lbl = 'natural'
        if simname == 'nref11n_v2_selfshield_z15': lbl = 'natural_v2'
        elif simname == 'nref11n_nref10f_selfshield_z6': lbl = 'forced'






        minr = 20
        maxr = 750
        #gd = where((array(R_c) > maxr)) [0]
        #ax1.plot(array(tt_arr)[gd], zeros(len(gd)), 'o', alpha = 0.4, color = clrs[sm], markersize = 5)
        #ax2.plot(array(tt_arr)[gd], zeros(len(gd)), 'o', alpha = 0.4, color = clrs[sm], markersize = 5)
        #ax3.plot(array(tt_arr)[gd], zeros(len(gd)), 'o', alpha = 0.4, color = clrs[sm], markersize = 5)
        #ax4.plot(array(tt_arr)[gd], zeros(len(gd)), 'o', alpha = 0.4, color = clrs[sm], markersize = 5)

        gd = where(array(R_c) < maxr) [0]
        ax1.plot(array(tt_arr)[gd], array(sfr_arr)[gd], '-', color = clrs[sm], label = lbl, linewidth = 3)
        ax2.plot(array(tt_arr)[gd], array(R_90)[gd], '-', color = clrs[sm], label = lbl,linewidth = 3)
        ax3.plot(array(tt_arr)[gd], array(sms_arr)[gd], '-', color = clrs[sm],label = lbl, linewidth = 3)
        ax4.plot(array(tt_arr)[gd], array(gms_arr)[gd], '-', color = clrs[sm],label = lbl, linewidth = 3)


        gd = where((array(R_c) < 15)) [0]
        if sm == 0: 
            ymn = 0
            ymx = 0.2
        if sm == 1:
            ymn = 0.8
            ymx = 1.0

        for g in gd:
            ax1.axvline(x = array(tt_arr)[g], ymin = ymn, ymax = ymx, alpha = 0.3, linewidth = 2, color = clrs[sm], zorder = 1)
            ax2.axvline(x = array(tt_arr)[g], ymin = ymn, ymax = ymx, alpha = 0.3, linewidth = 2, color = clrs[sm], zorder = 1)
            ax3.axvline(x = array(tt_arr)[g], ymin = ymn, ymax = ymx, alpha = 0.3, linewidth = 2, color = clrs[sm], zorder = 1)
            ax4.axvline(x = array(tt_arr)[g], ymin = ymn, ymax = ymx, alpha = 0.3, linewidth = 2, color = clrs[sm], zorder = 1)


        gd = where((array(R_c) > 141/0.7)) [0]
        if len(gd) > 0:

            if len(gd) == len(R_c):
                ax1.axvline(x = array(tt_arr)[gd[-1]],ymin = 0, ymax = 1,  alpha = 1.0, linestyle = 'dashed', linewidth = 2, color = 'black', zorder = 1)
                ax2.axvline(x = array(tt_arr)[gd[-1]],ymin = 0, ymax = 1,  alpha = 1.0, linestyle = 'dashed', linewidth = 2, color = 'black', zorder = 1)
                ax3.axvline(x = array(tt_arr)[gd[-1]],ymin = 0, ymax = 1,  alpha = 1.0, linestyle = 'dashed', linewidth = 2, color = 'black', zorder = 1)
                ax4.axvline(x = array(tt_arr)[gd[-1]],ymin = 0, ymax = 1,  alpha = 1.0, linestyle = 'dashed', linewidth = 2, color = 'black', zorder = 1)
            else:
                ax1.axvline(x = array(tt_arr)[gd[-1]],ymin = 0, ymax = 1,  alpha = 1.0, linestyle = 'dashed', linewidth = 2, color = clrs[sm], zorder = 1)
                ax2.axvline(x = array(tt_arr)[gd[-1]],ymin = 0, ymax = 1,  alpha = 1.0, linestyle = 'dashed', linewidth = 2, color = clrs[sm], zorder = 1)
                ax3.axvline(x = array(tt_arr)[gd[-1]],ymin = 0, ymax = 1,  alpha = 1.0, linestyle = 'dashed', linewidth = 2, color = clrs[sm], zorder = 1)
                ax4.axvline(x = array(tt_arr)[gd[-1]],ymin = 0, ymax = 1,  alpha = 1.0, linestyle = 'dashed', linewidth = 2, color = clrs[sm], zorder = 1)



        #ax1.set_ylim(0, 0.8)
    for ax in [ax1, ax2, ax3, ax4]:
        ax.legend(fontsize = 20, loc = 1)



    #ax.set_yscale('log')
    #ax2.set_ylim(0, 3)

    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xlabel('time (Gyr)', fontsize = 18)

    ax1.set_ylabel(r'star-formation rate (M$_{\odot}$ yr$^{-1}$)', fontsize = 18)
    ax2.set_ylabel(r'r$_{*,90}$ (kpc)', fontsize = 18)
    ax3.set_ylabel(r'M$_{*}$ (M$_{\odot}$)', fontsize = 18)
    ax4.set_ylabel(r'M$_{g}$ (M$_{\odot}$)', fontsize = 18)

    figname = fig_dir + '/SFR_n%.2i_f%.2i.png'%(sat_ns[0], sat_ns[1])
    fig1.tight_layout()
    fig1.savefig(figname, dpi = 500)

    figname = fig_dir + '/RE_n%.2i_f%.2i.png'%(sat_ns[0], sat_ns[1])
    fig2.tight_layout()
    fig2.savefig(figname, dpi = 500)


    figname = fig_dir + '/SM_n%.2i_f%.2i.png'%(sat_ns[0], sat_ns[1])
    fig3.tight_layout()
    fig3.savefig(figname, dpi = 500)

    figname = fig_dir + '/GM_n%.2i_f%.2i.png'%(sat_ns[0], sat_ns[1])
    fig4.tight_layout()
    fig4.savefig(figname, dpi = 500)


    print 'saved %s'%figname










