import astropy
from astropy.io import fits
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.major.size'] = 6

plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.major.size'] = 6



plt.ioff()
plt.close('all')



fig_dir = '/Users/rsimons/Dropbox/rcs_foggie/figures'
fits_dir = '/Users/rsimons/Dropbox/rcs_foggie/satellite_masses/all'







for s, simname in enumerate(['natural', 'nref11n_nref10f_selfshield_z6']):#, 'nref11n_selfshield_z15']):
    DDmin = 44
    DDmax = 500
    #if s == 1 : DDmax = 1050
    #else:       DDmax = 900
    #if s == 0: stn = 4
    #elif s == 1: stn = 3
    if simname == 'natural': nsats = 12
    else: nsats = 11
    for sat_n in arange(nsats):
        if sat_n !=9:

            fig2, (axes2, axes3, axes4) = plt.subplots(1,3, figsize = (21, 8))

            fig, axes = plt.subplots(1,1, figsize = (8, 8))

            dd_time = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/DD_time.npy')[()]




            zs, sm1, sm2, sm3, sm4 = [], [], [], [], []
            zs, ysm1, ysm2, ysm3, ysm4 = [], [], [], [], []
            zs, gm1, gm2, gm3, gm4 = [], [], [], [], []

            for d, DD in enumerate(arange(DDmin+1, DDmax, 1)):
                z = dd_time[1][where(dd_time[0] == DD)]

                data_p_fl = fits_dir + '/%s_DD%.4i_mass_sat%.2i.fits'%(simname, DD-1, sat_n)
                data_fl = fits_dir + '/%s_DD%.4i_mass_sat%.2i.fits'%(simname, DD, sat_n)

                if os.path.isfile(data_p_fl) & os.path.isfile(data_fl):
                    data_p = fits.open(data_p_fl)
                    data   = fits.open(data_fl)


                    r_cen = mean((data['distance'].data[1:], data['distance'].data[:-1]), axis = 0)
                    V_she = 4/3. * pi * (data['distance'].data[1:]**3. - data['distance'].data[:-1]**3.)

                    sM_den = (data['STARS_MASS'].data[1:] - data['STARS_MASS'].data[:-1])/V_she
                    sM_den_p = (data_p['STARS_MASS'].data[1:] - data_p['STARS_MASS'].data[:-1])/V_she

                    ysM_den = (data['STARS_YOUNGMASS'].data[1:] - data['STARS_YOUNGMASS'].data[:-1])/V_she
                    ysM_den_p = (data_p['STARS_YOUNGMASS'].data[1:] - data_p['STARS_YOUNGMASS'].data[:-1])/V_she


                    gM_den = (data['GAS_TOT'].data[1:] - data['GAS_TOT'].data[:-1])/V_she
                    gM_den_p = (data_p['GAS_TOT'].data[1:] - data_p['GAS_TOT'].data[:-1])/V_she
                    
                    if d == 0:
                        axes.plot(r_cen, sM_den, color = 'red')
                        axes.plot(r_cen, gM_den, color = 'blue')


                    zs.append(z)

                    sm1.append(sM_den[0])
                    sm2.append(sM_den[2])
                    sm3.append(sM_den[11])
                    sm4.append(sM_den[15])

                    ysm1.append(ysM_den[0])
                    ysm2.append(ysM_den[2])
                    ysm3.append(ysM_den[11])
                    ysm4.append(ysM_den[15])

                    gm1.append(gM_den[0])
                    gm2.append(gM_den[2])
                    gm3.append(gM_den[11])
                    gm4.append(gM_den[15])







            fs = 25
            #if s == 0:
            axes2.annotate('R = 0.4 kpc', (5.8, 1.7 * sm1[0]), color = 'black',  fontweight = 'bold', fontsize = fs, zorder = 10)
            axes2.annotate('1 kpc',       (5.8, 1.7 * sm2[0]), color = 'black',  fontweight = 'bold', fontsize = fs, zorder = 10)
            axes2.annotate('3 kpc',       (5.8, 1.7 * sm3[0]), color = 'black',  fontweight = 'bold', fontsize = fs, zorder = 10)


            #if s == 0:
            clr_s = ['darkred', 'darkred', 'darkred']
            clr_g = ['darkblue', 'darkblue', 'darkblue']
            #else:
            #    clr_s = ['red', 'red', 'red']
            #    clr_g = ['blue', 'blue', 'blue']

            alps = [1.0, 0.7, 0.5]

            axes2.plot(array(zs), array(sm1), linewidth = 4,color = clr_s[0], alpha = alps[0], linestyle = '-' , zorder = 9)
            axes2.plot(array(zs), array(sm2), linewidth = 4,color = clr_s[1], alpha = alps[1], linestyle = '--', zorder = 9)
            axes2.plot(array(zs), array(sm3), linewidth = 4,color = clr_s[2], alpha = alps[2], linestyle = ':' , zorder = 9)

            axes3.plot(array(zs), array(gm1), linewidth = 4, color = clr_g[0],  alpha = alps[0], linestyle = '-' , zorder = 9)
            axes3.plot(array(zs), array(gm2), linewidth = 4, color = clr_g[1],  alpha = alps[1], linestyle = '--', zorder = 9)
            axes3.plot(array(zs), array(gm3), linewidth = 4, color = clr_g[2],  alpha = alps[2], linestyle = ':' , zorder = 9)

            axes4.plot(array(zs), array(gm1)/(array(gm1) + sm1), linewidth = 4, color = 'black',  alpha = alps[0], linestyle = '-' , zorder = 9)
            axes4.plot(array(zs), array(gm2)/(array(gm2) + sm2), linewidth = 4, color = 'black',  alpha = alps[1], linestyle = '--', zorder = 9)
            axes4.plot(array(zs), array(gm3)/(array(gm3) + sm3), linewidth = 4, color = 'black',  alpha = alps[2], linestyle = ':' , zorder = 9)


            #axes2.plot(array(zs), array(sm4), linewidth = 4,color = 'red', alpha = 0.3, linestyle = '--')
            #axes3.plot(array(zs), array(gm4), linewidth = 4,color = 'blue', alpha = 0.3, linestyle = '--')
            #axes2.annotate('5 kpc', (1.52, 1.5 * sm4[0]),color = 'black',   fontweight = 'bold', fontsize = fs,)


            #if s == 0:
            axes3.annotate('Satellite %s'%sat_n, (0.5, 0.92), ha = 'center', xycoords = 'figure fraction', fontsize = 40, fontweight = 'bold')


            for ax in [axes2, axes3, axes4]:
                ax.set_xlim(6, 1.4)
                ax.set_xlabel('redshift', fontsize = 25)




            for ax in [axes, axes2, axes3]:
                ax.set_yscale('log')



            axes.set_ylabel(r'$\rho$ (M$_{\odot}$ kpc$^{-3}$)', fontsize = 25)
            axes2.set_ylabel(r'$\rho_*$ (M$_{\odot}$ kpc$^{-3}$)', fontsize = 25)
            axes3.set_ylabel(r'$\rho_g$ (M$_{\odot}$ kpc$^{-3}$)', fontsize = 25)

            axes4.set_ylabel(r'gas fraction', fontsize = 25)





            axes4.set_ylim(0,1)


            ylm_mn = min((axes2.get_ylim()[0], axes3.get_ylim()[0]))
            ylm_mx = max((axes2.get_ylim()[1], axes3.get_ylim()[1]))

            #ylm_mn = 1.e2
            #ylm_mx = 1.e8

            for ax in [axes2, axes3]:
                ax.set_ylim(ylm_mn, ylm_mx)



            #axes.set_xscale('log')


            #axes.set_xscale('linear')
            axes.set_xlabel('distance from center (kpc)', fontsize = 25)


            axes.set_xticks(arange(0, 20), minor = True)
            axes.set_xticks([0, 5, 10, 15])
            axes.set_xticklabels(['0', '5', '10', '15'])

            #axes1[1].set_yscale('symlog')


            #axes1[1].set_ylabel(r'$\dot{\rho_*}$ (M$_{\odot}$ kpc$^{-3}$ yr$^{-1}$)', fontsize = 20)


            #for ax in [axes1[1], axes2[1]]:
            #    axes1[1].axhline(y = 0, color = 'black', linestyle = '--')

            #for fig in [fig1, fig2]:
            axes.annotate('gas', (0.7, 0.9), xycoords = 'axes fraction', color = 'blue', fontweight = 'bold', fontsize = 40)
            axes.annotate('stars', (0.7, 0.8), xycoords = 'axes fraction', color = 'red', fontweight = 'bold', fontsize = 40)

            fig.subplots_adjust(left = 0.12, right = 0.95, wspace = 0.30)

            fig.savefig(fig_dir + '/profiles/%s_masses_sat%.2i.png'%(simname, sat_n), dpi = 300)
            fig2.subplots_adjust(left = 0.05, right = 0.98, hspace = 0.33)
            fig2.savefig(fig_dir + '/time_profiles/%s_masses_time_sat%.2i.png'%(simname, sat_n), dpi = 300)





            # Make movie frames
            if False:
                for d, DD in enumerate(arange(DDmin+1, DDmax, 1)):
                    figxi, (axes2xi, axes3xi) = plt.subplots(1,2, figsize = (15, 8))

                    alps = [1.0, 0.7, 0.5]
                    axes2xi.plot(array(zs), array(sm1), linewidth = 4,color = clr_s[0], alpha = alps[0], linestyle = '-' , zorder = 9)
                    axes2xi.plot(array(zs), array(sm2), linewidth = 4,color = clr_s[1], alpha = alps[1], linestyle = '--', zorder = 9)
                    axes2xi.plot(array(zs), array(sm3), linewidth = 4,color = clr_s[2], alpha = alps[2], linestyle = ':' , zorder = 9)

                    axes3xi.plot(array(zs), array(gm1), linewidth = 4, color = clr_g[0],  alpha = alps[0], linestyle = '-' , zorder = 9)
                    axes3xi.plot(array(zs), array(gm2), linewidth = 4, color = clr_g[1],  alpha = alps[1], linestyle = '--', zorder = 9)
                    axes3xi.plot(array(zs), array(gm3), linewidth = 4, color = clr_g[2],  alpha = alps[2], linestyle = ':' , zorder = 9)


                    #axes3.plot(array(zs[::-1]),  array(gm1[0]) - cumsum(array(ysm1)[::-1]), linewidth = 4, color = 'darkgreen', alpha = 0.4, linestyle = '--')
                    #axes3.plot(array(zs[::-1]),  array(gm2[0]) - cumsum(array(ysm2)[::-1]), linewidth = 4, color = 'darkblue',  alpha = 0.4, linestyle = '--')
                    #axes3.plot(array(zs[::-1]),  array(gm3[0]) - cumsum(array(ysm3)[::-1]), linewidth = 4, color = 'darkred',   alpha = 0.4, linestyle = '--')

                    fs = 25
                    axes2xi.annotate('R = 0.4 kpc', (1.52, 1.7 * sm1[0]), color = 'black',  fontweight = 'bold', fontsize = fs, zorder = 10)
                    axes2xi.annotate('1 kpc',       (1.52, 1.7 * sm2[0]), color = 'black',  fontweight = 'bold', fontsize = fs, zorder = 10)
                    axes2xi.annotate('3 kpc',       (1.52, 1.7 * sm3[0]), color = 'black',  fontweight = 'bold', fontsize = fs, zorder = 10)
                    z = dd_time[1][where(dd_time[0] == DD)]
                    axes2xi.axvline(x = z, color = 'black', linewidth = 3, zorder = 1)
                    axes3xi.axvline(x = z, color = 'black', linewidth = 3, zorder = 1)



                    axes3xi.annotate('Satellite %s'%sat_n, (0.5, 0.92), ha = 'center', xycoords = 'figure fraction', fontsize = 40, fontweight = 'bold')


                    for ax in [axes2xi, axes3xi]:
                        ax.set_xlim(1.55, 0.9)
                        ax.set_xlabel('redshift', fontsize = 25)
                        ax.set_yscale('log')

                    axes2xi.set_ylabel(r'$\rho_*$ (M$_{\odot}$ kpc$^{-3}$)', fontsize = 25)
                    axes3xi.set_ylabel(r'$\rho_g$ (M$_{\odot}$ kpc$^{-3}$)', fontsize = 25)




                    ylm_mn = min((axes2xi.get_ylim()[0], axes3xi.get_ylim()[0]))
                    ylm_mx = max((axes2xi.get_ylim()[1], axes3xi.get_ylim()[1]))


                    for ax in [axes2xi, axes3xi]:
                        ax.set_ylim(ylm_mn, ylm_mx)


                    figxi.subplots_adjust(left = 0.08, right = 0.98, wspace = 0.30)
                    figxi.savefig(fig_dir + '/masses_time_DD/%s_masses_time_sat%.2i_DD%.4i_new.png'%(simname, sat_n, DD), dpi = 300)
















            plt.close('all')
                
    #fig2.subplots_adjust(left = 0.05, right = 0.98, hspace = 0.33)
    #fig2.savefig(fig_dir + '/time_profiles/both_masses_time_F4N3.png', dpi = 300)
    #fig2.savefig(fig_dir + '/time_profiles/both_masses_time_F4N3.png', dpi = 300)




    
    
    
    
    