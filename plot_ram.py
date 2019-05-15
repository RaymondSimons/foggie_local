import astropy
from astropy.io import fits
plt.ioff()
plt.close('all')



fig1, axes1 = plt.subplots(2,3, figsize = (12, 8))
fig2, axes2 = plt.subplots(2,3, figsize = (12, 8))
fig3, axes3 = plt.subplots(2,3, figsize = (12, 8))


fig1_2, axes1_2 = plt.subplots(2,3, figsize = (12, 8))
fig2_2, axes2_2 = plt.subplots(2,3, figsize = (12, 8))
fig3_2, axes3_2 = plt.subplots(2,3, figsize = (12, 8))


fig1_3, axes1_3 = plt.subplots(2,3, figsize = (12, 8))
fig2_3, axes2_3 = plt.subplots(2,3, figsize = (12, 8))
fig3_3, axes3_3 = plt.subplots(2,3, figsize = (12, 8))


simnames =  ['natural',
        'natural_v2',
        'natural_v3',
        'natural_v4',
        'nref11c_nref9f',
        'nref11n_nref10f']


clrs = ['blue', 'navy', 'darkblue', 'royalblue', 'green', 'red']

'''
simnames =  ['natural',
            'nref11c_nref9f',
            'nref11n_nref10f']

simnames =  ['natural',
            'nref11n_nref10f']

clrs = ['blue', 'red', 'green']
'''
DD_to_t = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/DD_time.npy')[()]


ls = ['-', ':', '--','-.', '-', '-']
lw = [2,2,2,2,2, 2]
alp = [0.8, 0.8, 0.8, 0.8, 0.8,  0.8]

tmin = [1.3, 1.4, 1.3, 1.3, 1.2, 1.3]
tmax = [5.65, 3.5, 2.5, 2.2, 1.6, 1.9]

axs = 'y'
for sat in arange(6):
    ax1 = axes1.ravel()[sat]
    ax2 = axes2.ravel()[sat]
    ax3 = axes3.ravel()[sat]


    ax1_2 = axes1_2.ravel()[sat]
    ax2_2 = axes2_2.ravel()[sat]
    ax3_2 = axes3_2.ravel()[sat]


    ax1_3 = axes1_3.ravel()[sat]
    ax2_3 = axes2_3.ravel()[sat]
    ax3_3 = axes3_3.ravel()[sat]

    for s, simname in enumerate(simnames):
        data = fits.open('/Users/rsimons/Dropbox/rcs_foggie/ram_pressure/percentiles/%s_ram_percentiles.fits'%simname)
        DD = data['DD'].data['DD']
        td = array([DD_to_t[2][where(DD_to_t[0] == d)[0][0]] for d in DD])
        data = data['SAT_%.2i'%sat]

        gd = where((~isnan(data.data['density_02'])) & (td < tmax[sat]) & (td > tmin[sat]))[0]

        '''
        lws = [1, 1,2,1,1 ]
        ls = ['--', '-', '-', '-', '--']
        lalp = [0.2, 0.3, 1.0, 0.3, 0.2]
        for p, perc in enumerate(['02', '16', '50', '84', '98']):
            ax1.plot(td[gd], data.data['density_%s'%perc][gd], color = clrs[s], linewidth = lws[p], linestyle = ls[p], alpha = lalp[p])
            ax2.plot(td[gd], data.data['x_velocity_%s'%perc][gd], color = clrs[s],  linewidth = lws[p], linestyle = ls[p], alpha = lalp[p])
            ax3.plot(td[gd], data.data['x_ram_%s'%perc][gd], color = clrs[s],  linewidth = lws[p], linestyle = ls[p], alpha = lalp[p])
        '''
        if (simname == 'natural') | (simname == 'nref11n_nref10f'):
            for p, perc in enumerate(['02', '16', '50', '84', '98']):
                ax1.fill_between(td[gd], data.data['density_98'][gd],     data.data['density_02'][gd],     color = clrs[s], alpha = 0.3)#, linewidth = lws[p], linestyle = ls[p], alpha = lalp[p])
                ax2.fill_between(td[gd], data.data['%s_velocity_98'%axs][gd],  data.data['%s_velocity_02'%axs][gd],  color = clrs[s], alpha = 0.3)#,  linewidth = lws[p], linestyle = ls[p], alpha = lalp[p])
                ax3.fill_between(td[gd], data.data['%s_ram_98'%axs][gd],       data.data['%s_ram_02'%axs][gd],       color = clrs[s], alpha = 0.3)#,  linewidth = lws[p], linestyle = ls[p], alpha = lalp[p])
     

            lws = [2]
            lss = ['-']
            lalp = [1.0]
            for p, perc in enumerate(['50']):
                if sat == 0:
                    ax1.plot(td[gd], data.data['density_%s'%perc][gd], color = clrs[s], linewidth = lws[p], linestyle = lss[p], alpha = lalp[p], label = simname.replace('_', '-'))
                    ax2.plot(td[gd], data.data['%s_velocity_%s'%(axs, perc)][gd], color = clrs[s],  linewidth = lws[p], linestyle = lss[p], alpha = lalp[p], label = simname.replace('_', '-'))
                    ax3.plot(td[gd], data.data['%s_ram_%s'%(axs, perc)][gd], color = clrs[s],  linewidth = lws[p], linestyle = lss[p], alpha = lalp[p], label = simname.replace('_', '-'))
                    for ax in [ax1, ax2, ax3]:
                        ax.legend(loc = 4)

                else:
                    ax1.plot(td[gd], data.data['density_%s'%perc][gd], color = clrs[s], linewidth = lws[p], linestyle = ls[p], alpha = lalp[p])
                    ax2.plot(td[gd], data.data['%s_velocity_%s'%(axs, perc)][gd], color = clrs[s],  linewidth = lws[p], linestyle = ls[p], alpha = lalp[p])
                    ax3.plot(td[gd], data.data['%s_ram_%s'%(axs, perc)][gd], color = clrs[s],  linewidth = lws[p], linestyle = ls[p], alpha = lalp[p])

        if sat ==0:
            ax1_2.plot(td[gd], (data.data['density_84'][gd]     - data.data['density_16'][gd]    )/2., color = clrs[s], linewidth  = 1.0,linestyle =  '-', alpha = 1.0,label = simname.replace('_', '-'))
            ax2_2.plot(td[gd], (data.data['%s_velocity_84'%axs][gd]  - data.data['%s_velocity_16'%axs][gd] )/2., color = clrs[s],  linewidth = 1.0, linestyle = '-', alpha = 1.0,label = simname.replace('_', '-'))
            ax3_2.plot(td[gd], (data.data['%s_ram_84'%axs][gd]       - data.data['%s_ram_16'%axs][gd]      )/2., color = clrs[s],  linewidth = 1.0, linestyle = '-', alpha = 1.0,label = simname.replace('_', '-'))
            ax1_3.plot(td[gd], (data.data['density_98'][gd]     - data.data['density_02'][gd]    )/2., color = clrs[s], linewidth  = 1.0,linestyle =  '-', alpha = 1.0,label = simname.replace('_', '-'))
            ax2_3.plot(td[gd], (data.data['%s_velocity_98'%axs][gd]  - data.data['%s_velocity_02'%axs][gd] )/2., color = clrs[s],  linewidth = 1.0, linestyle = '-', alpha = 1.0,label = simname.replace('_', '-'))
            ax3_3.plot(td[gd], (data.data['%s_ram_98'%axs][gd]       - data.data['%s_ram_02'%axs][gd]      )/2., color = clrs[s],  linewidth = 1.0, linestyle = '-', alpha = 1.0,label = simname.replace('_', '-'))
            for ax in [ax1_2, ax2_2, ax3_2, ax1_3, ax2_3, ax3_3]:
                ax.legend(loc = 1)

        else:
            ax1_2.plot(td[gd], (data.data['density_84'][gd]     - data.data['density_16'][gd]    )/2., color = clrs[s], linewidth  = 1.0,linestyle = '-', alpha = 1.0)
            ax2_2.plot(td[gd], (data.data['%s_velocity_84'%axs][gd]  - data.data['%s_velocity_16'%axs][gd] )/2., color = clrs[s],  linewidth = 1.0, linestyle ='-', alpha = 1.0)
            ax3_2.plot(td[gd], (data.data['%s_ram_84'%axs][gd]       - data.data['%s_ram_16'%axs][gd]      )/2., color = clrs[s],  linewidth = 1.0, linestyle ='-', alpha = 1.0)
            ax1_3.plot(td[gd], (data.data['density_98'][gd]     - data.data['density_02'][gd]    )/2., color = clrs[s], linewidth  = 1.0,linestyle = '-', alpha = 1.0)
            ax2_3.plot(td[gd], (data.data['%s_velocity_98'%axs][gd]  - data.data['%s_velocity_02'%axs][gd] )/2., color = clrs[s],  linewidth = 1.0, linestyle ='-', alpha = 1.0)
            ax3_3.plot(td[gd], (data.data['%s_ram_98'%axs][gd]       - data.data['%s_ram_02'%axs][gd]      )/2., color = clrs[s],  linewidth = 1.0, linestyle ='-', alpha = 1.0)


        #ax.fill_between(DD[gd], data.data['density_02'][gd], data.data['density_98'][gd], color = clrs[s], alpha = 0.3)
        #ax2.fill_between(DD[gd], data.data['x_velocity_02'][gd], data.data['x_velocity_98'][gd], color = clrs[s], alpha = 0.3)

    ax1.set_yscale('log')
    ax3.set_yscale('log')

    fs = 14

    for ax in [ax1, ax2, ax3, ax1_2, ax3_2, ax3_2, ax1_3, ax2_3, ax3_3]:
        if sat > 2:
            ax.set_xlabel('time (Gyr)', fontsize = fs)
        ax.annotate('Satellite %i'%sat, (0.05, 0.85), xycoords = 'axes fraction', fontsize = 16)

    if sat%3 == 0:
        ax1.set_ylabel('gas density (g cm$^{-3}$)', fontsize = fs)
        ax2.set_ylabel('gas velocity (km s$^{-1}$)', fontsize = fs)
        ax3.set_ylabel('Ram Pressure (dyn cm$^{-2}$)', fontsize = fs)

        ax1_2.set_ylabel('$\sigma_{68}$ gas density (g cm$^{-3}$)', fontsize = fs)
        ax2_2.set_ylabel('$\sigma_{68}$ gas velocity (km s$^{-1}$)', fontsize = fs)
        ax3_2.set_ylabel('$\sigma_{68}$ Ram Pressure (dyn cm$^{-2}$)', fontsize = fs)

        ax1_3.set_ylabel('$\sigma_{95}$ gas density (g cm$^{-3}$)', fontsize = fs)
        ax2_3.set_ylabel('$\sigma_{95}$ gas velocity (km s$^{-1}$)', fontsize = fs)
        ax3_3.set_ylabel('$\sigma_{95}$ Ram Pressure (dyn cm$^{-2}$)', fontsize = fs)




fig1.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/rams/dens.png', dpi = 300)
fig2.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/rams/vel.png', dpi = 300)
fig3.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/rams/ram.png', dpi = 300)





fig1_2.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/rams/all/dens_68.png', dpi = 300)
fig2_2.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/rams/all/vel_68.png', dpi = 300)
fig3_2.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/rams/all/ram_68.png', dpi = 300)



fig1_3.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/rams/all/dens_95.png', dpi = 300)
fig2_3.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/rams/all/vel_95.png', dpi = 300)
fig3_3.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/rams/all/ram_95.png', dpi = 300)







