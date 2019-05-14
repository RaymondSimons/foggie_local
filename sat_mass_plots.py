import glob
from glob import glob
from astropy.io import fits
from scipy import interpolate
plt.rcParams['text.usetex'] = True
plt.ioff()
plt.close('all')

DD_to_t = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/DD_time.npy')[()]


DDs_use = arange(49, 1000)

clrs = ['blue', 'navy', 'darkblue', 'royalblue', 'red', 'green']
ls = ['-', ':', '--','-.', '-', '-']
lw = [2,2,2,2,2, 2]
alp = [0.8, 0.8, 0.8, 0.8, 0.8,  0.8]
#for s, sim in enumerate(array(['forced', 'natural(v1)', 'natural(v2)', 'natural(v3)', 'natural(v4)'])):
for sat_n in arange(0,6):
    fig, axes = plt.subplots(2,2, figsize = (10, 10))
    fig2, ax2 = plt.subplots(1,1, figsize = (12, 6))
    sims = array(['natural(v1)', 'natural(v2)', 'natural(v3)', 'natural(v4)', 'nref11c-nref9f','nref11n-nref10f'])
    ts = nan * zeros((len(sims), len(DDs_use)))
    ms = nan * zeros((len(sims), len(DDs_use)))
    mg = nan * zeros((len(sims), len(DDs_use)))
    sf = nan * zeros((len(sims), len(DDs_use)))
    dm = nan * zeros((len(sims), len(DDs_use)))



    for s, sim in enumerate(sims):
        print sim
        for d, DD in enumerate(DDs_use):
            if sim == 'natural(v1)': 
                simname = 'natural'
                simfname = 'natural'

            if sim == 'natural(v2)': 
                simname = 'natural_v2'
                simfname = 'nref11n_v2_selfshield_z15'

            if sim == 'natural(v3)': 
                simname = 'natural_v3'
                simfname = 'nref11n_v3_selfshield_z15'

            if sim == 'natural(v4)': 
                simname = 'natural_v4'
                simfname = 'nref11n_v4_selfshield_z15'

            if sim == 'nref11c-nref9f': 
                simname = 'nref11c_nref9f'
                simfname = 'nref11c_nref9f'

            if sim == 'nref11n-nref10f': 
                simname = 'nref11n_nref10f'
                simfname = 'nref11n_nref10f'

            #fl = glob('/Users/rsimons/Dropbox/rcs_foggie/cenmass/%s_DD%.4i_mass.fits'%(simname, DD))
            fl = glob('/Users/rsimons/Dropbox/rcs_foggie/satellite_masses/%s/%s_DD%.4i_mass.fits'%(simname, simfname, DD))

            if len(fl) > 0:
                a = fits.open(fl[0])
                a = a['SAT_%.2i'%sat_n].data
                gd = where(DD_to_t[0] == DD)
                t = DD_to_t[2][gd][0]
                frac_mass = a['STARS_MASS']/a['STARS_MASS'][-1]
                ts[s, d] = t
                ms[s, d] = a['STARS_MASS'][-1]
                mg[s, d] = a['GAS_TOT'][-1]
                sf[s, d] = a['STARS_YOUNGMASS'][-1]/2.e7
                dm[s, d] = a['DARK_MATTER'][-1]
                #f = interpolate.interp1d(frac_mass, a['RADIUS'].data)
                #R_90.append(f(0.9))




        axes[0,0].plot(ts[s,:], ms[s,:], label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        axes[0,1].plot(ts[s,:], dm[s,:], label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        axes[1,0].plot(ts[s,:], mg[s,:], label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        #axes[1,1].plot(ts, R_90, label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        ax2.plot(ts[s,:], sf[s,:], label = sim, color = clrs[s], linestyle = ls[s], linewidth = 1, alpha = alp[s])



    axes[0,0].fill_between(ts[s,:], np.min(ms[0:4,:], axis = 0), np.max(ms[0:4,:], axis = 0), color = 'blue', alpha = 0.2)
    axes[0,1].fill_between(ts[s,:], np.min(dm[0:4,:], axis = 0), np.max(dm[0:4,:], axis = 0), color = 'blue', alpha = 0.2)
    axes[1,0].fill_between(ts[s,:], np.min(mg[0:4,:], axis = 0), np.max(mg[0:4,:], axis = 0), color = 'blue', alpha = 0.2)

    axes[0,0].legend(loc = 2)
    ax2.legend(loc = 2)


    fs = 12

    axes[0,0].set_ylabel('M$_*$ (M$_{\odot}$)', fontsize = fs)
    axes[0,1].set_ylabel('M$_{DM}$ (M$_{\odot}$)', fontsize = fs)
    axes[1,0].set_ylabel('M$_{g}$ (M$_{\odot}$)', fontsize = fs)
    axes[1,1].set_ylabel(r'r$_{*,90}$ (kpc)', fontsize = fs)
    ax2.set_ylabel('star formation rate (M$_{\odot}$ yr$^{-1}$)', fontsize = fs)




    for ax in axes.ravel(): ax.set_xlabel('time (Gyr)', fontsize = fs)
    ax2.set_xlabel('time (Gyr)', fontsize = fs)
    fig.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/butterfly/%i_mass.png'%sat_n, dpi = 300)
    fig2.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/butterfly/%i_SFR.png'%sat_n, dpi = 300)

    
    






