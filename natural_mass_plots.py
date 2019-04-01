import glob
from glob import glob
from astropy.io import fits
plt.rcParams['text.usetex'] = True
plt.ioff()
plt.close('all')

DD_to_t = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/DD_time.npy')[()]


DDs = arange(44, 1000)
fig, axes = plt.subplots(2,2, figsize = (10, 10))
clrs = ['red', 'blue', 'green']
for s, sim in enumerate(array(['v2', 'v3', 'v4'])):
    print sim
    ts = []
    ms = []
    mg = []
    sf = []
    dm = []
    for DD in DDs:
        fl = glob('/Users/rsimons/Dropbox/rcs_foggie/cenmass/nref11n_%s_selfshield_z15_DD%.4i_mass.fits'%(sim, DD))
        if len(fl) > 0:
            a = fits.open(fl[0])
            gd = where(DD_to_t[0] == DD)
            t = DD_to_t[2][gd][0]
            ts.append(t)
            ms.append(a['STARS_MASS'].data[-1])
            mg.append(a['GAS_TOT'].data[-1])
            sf.append(a['STARS_YOUNGMASS'].data[-1]/2.e7)
            dm.append(a['DARK_MATTER'].data[-1])

    axes[0,0].plot(ts, ms, label = sim, color = clrs[s])
    axes[0,1].plot(ts, dm, label = sim, color = clrs[s])
    axes[1,0].plot(ts, mg, label = sim, color = clrs[s])
    axes[1,1].plot(ts, sf, label = sim, color = clrs[s])

axes[0,0].legend(loc = 2)


fs = 12

axes[0,0].set_ylabel('M$_*$ (M$_{\odot}$)', fontsize = fs)
axes[0,1].set_ylabel('M$_{DM}$ (M$_{\odot}$)', fontsize = fs)
axes[1,0].set_ylabel('M$_{g}$ (M$_{\odot}$)', fontsize = fs)
axes[1,1].set_ylabel('star-formation rate (M$_{\odot}$ yr$^{-1}$)', fontsize = fs)



for ax in axes.ravel(): ax.set_xlabel('time (Gyr)', fontsize = fs)
fig.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/naturals_mass.png', dpi = 300)

