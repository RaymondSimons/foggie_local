import glob
from glob import glob
from astropy.io import fits
from scipy import interpolate
plt.rcParams['text.usetex'] = True
plt.ioff()
plt.close('all')

DD_to_t = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/DD_time.npy')[()]


#DDs = arange(44, 800)
DDs = arange(44, 800)
sat_ns = arange(0, 12)

clrs = ['blue', 'navy', 'darkblue', 'royalblue', 'red']

ls = ['-', '-', ':', '--','-.']
lw = [2,2,2,2,2]
alp = [0.8, 0.8, 0.8, 0.8, 0.8]
#for s, sim in enumerate(array(['forced', 'natural(v1)', 'natural(v2)', 'natural(v3)', 'natural(v4)'])):
for sat_n in sat_ns:
    print sat_n
    fig, axes = plt.subplots(2,2, figsize = (10, 10))
    fig2, ax2 = plt.subplots(1,1, figsize = (12, 6))
    for s, sim in enumerate(array(['natural(v1)', 'natural(v2)', 'natural(v3)', 'natural(v4)'])):
        ts = []
        ms = []
        mg = []
        sf = []
        dm = []
        R_90 = []
        #if sim == 'forced': DDs_use = arange(400, 1000)
        #elif sim == 'natural(v1)': DDs_use = arange(44, 800)
        #else: DDs_use = DDs
        DDs_use = DDs
        for DD in DDs_use:
            if sim == 'forced': simname = 'nref11n_nref10f'
            if sim == 'natural(v1)': 
                simname = 'natural'
                dirname = 'natural'
                ind_use = 50

            if sim == 'natural(v2)': 
                dirname = 'natural_v2'
                simname = 'nref11n_v2_selfshield_z15'
                ind_use = 0

            if sim == 'natural(v3)': 
                dirname = 'natural_v3'
                simname = 'nref11n_v3_selfshield_z15'
                ind_use = 0

            if sim == 'natural(v4)': 
                dirname = 'natural_v4'
                simname = 'nref11n_v4_selfshield_z15'
                ind_use = 0

            #fl = glob('/Users/rsimons/Dropbox/rcs_foggie/cenmass/%s_DD%.4i_mass.fits'%(simname, DD))
            fl_name = '/Users/rsimons/Dropbox/rcs_foggie/satmass/%s/%s_DD%.4i_mass_sat%.2i.fits'%(dirname, simname, DD, sat_n)
            fl = glob(fl_name)
            if len(fl) > 0:
                try:
                    a = fits.open(fl[0])
                    gd = where(DD_to_t[0] == DD)
                    t = DD_to_t[2][gd][0]
                    frac_mass = a['STARS_MASS'].data[:]/a['STARS_MASS'].data[-1]
                    ts.append(t)
                    ms.append(a['STARS_MASS'].data[ind_use])
                    mg.append(a['GAS_TOT'].data[ind_use])
                    sf.append(a['STARS_YOUNGMASS'].data[ind_use]/2.e7)
                    dm.append(a['DARK_MATTER'].data[ind_use])
                    #f = interpolate.interp1d(frac_mass, a['DISTANCE'].data)
                    #R_90.append(f(0.9))
                except:
                    pass



        axes[0,0].plot(ts, ms, label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        axes[0,1].plot(ts, dm, label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        axes[1,0].plot(ts, mg, label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        #axes[1,1].plot(ts, R_90, label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        ax2.plot(ts, sf, label = sim, color = clrs[s], linestyle = ls[s], linewidth = 1, alpha = alp[s])

    axes[0,0].legend(loc = 2)
    ax2.legend(loc = 2)


    fs = 12

    axes[0,0].set_ylabel('M$_*$ (M$_{\odot}$)', fontsize = fs)
    axes[0,1].set_ylabel('M$_{DM}$ (M$_{\odot}$)', fontsize = fs)
    axes[1,0].set_ylabel('M$_{g}$ (M$_{\odot}$)', fontsize = fs)
    axes[1,1].set_ylabel(r'r$_{*,90}$ (kpc)', fontsize = fs)
    ax2.set_ylabel('star formation rate (M$_{\odot}$ yr$^{-1}$)', fontsize = fs)


    for ax in axes.ravel(): ax.set_xlim(0, 5.5)
    ax2.set_xlim(0, 5.5)

    axes[1,1].axis('off')



    for ax in axes.ravel(): ax.set_xlabel('time (Gyr)', fontsize = fs)
    ax2.set_xlabel('time (Gyr)', fontsize = fs)
    fig.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/butterfly_sats/%i_mass.png'%sat_n, dpi = 300)
    fig2.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/butterfly_sats/%i_SFR.png'%sat_n, dpi = 300)

    plt.close('all')





