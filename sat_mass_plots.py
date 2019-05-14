import glob
from glob import glob
from astropy.io import fits
from scipy import interpolate
plt.rcParams['text.usetex'] = True
plt.ioff()
plt.close('all')

DD_to_t = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/DD_time.npy')[()]


DDs_use = arange(100, 1000)

clrs = ['blue', 'navy', 'darkblue', 'royalblue', 'red', 'green']
ls = ['-', ':', '--','-.', '-', '-']
lw = [2,2,2,2,2, 2]
alp = [0.8, 0.8, 0.8, 0.8, 0.8,  0.8]
#for s, sim in enumerate(array(['forced', 'natural(v1)', 'natural(v2)', 'natural(v3)', 'natural(v4)'])):


sats = arange(0,6)
sims = array(['natural(v1)', 'natural(v2)', 'natural(v3)', 'natural(v4)', 'nref11c-nref9f','nref11n-nref10f'])

ts = nan * zeros((len(sims), len(DDs_use), len(sats)))
ms = nan * zeros((len(sims), len(DDs_use), len(sats)))
mg = nan * zeros((len(sims), len(DDs_use), len(sats)))
sf = nan * zeros((len(sims), len(DDs_use), len(sats)))
dm = nan * zeros((len(sims), len(DDs_use), len(sats)))

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


        fl = glob('/Users/rsimons/Dropbox/rcs_foggie/satellite_masses/%s/%s_DD%.4i_mass.fits'%(simname, simfname, DD))

        if len(fl) > 0:
            data = fits.open(fl[0])
            for sat_n in sats:
                a = data['SAT_%.2i'%sat_n].data
                gd = where(DD_to_t[0] == DD)
                t = DD_to_t[2][gd][0]
                frac_mass = a['STARS_MASS']/a['STARS_MASS'][-1]
                ts[s, d, sat_n] = t
                ms[s, d, sat_n] = a['STARS_MASS'][-1]
                mg[s, d, sat_n] = a['GAS_TOT'][-1]
                sf[s, d, sat_n] = a['STARS_YOUNGMASS'][-1]/2.e7
                dm[s, d, sat_n] = a['DARK_MATTER'][-1]


for s, sim in enumerate(sims):
    master_hdulist = []
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)

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




    for sat_n in sats:
        cols = []
        cols.append(fits.Column(name = 'ts', array =  np.array(ts[s, :, sat_n]), format = 'D'))
        cols.append(fits.Column(name = 'ms', array =  np.array(ms[s, :, sat_n]), format = 'D'))
        cols.append(fits.Column(name = 'mg', array =  np.array(mg[s, :, sat_n]), format = 'D'))
        cols.append(fits.Column(name = 'sf', array =  np.array(sf[s, :, sat_n]), format = 'D'))
        cols.append(fits.Column(name = 'dm', array =  np.array(dm[s, :, sat_n]), format = 'D'))

        cols = fits.ColDefs(cols)


        master_hdulist.append(fits.BinTableHDU.from_columns(cols, name = 'SAT_%.2i'%sat_n))

    fits_name = '/Users/rsimons/Dropbox/rcs_foggie/satellite_masses/all/%s_all.fits'%(simname)
    thdulist = fits.HDUList(master_hdulist)
    print ('\tSaving to ' + fits_name)
    thdulist.writeto(fits_name, overwrite = True)



for sat_n in sats:
    fig, ax1 = plt.subplots(1,3, figsize = (15, 5))
    fig2, ax2 = plt.subplots(1,1, figsize = (12, 6))

    ms_bounds = nan*zeros((len(DDs_use), 2))
    dm_bounds = nan*zeros((len(DDs_use), 2))
    mg_bounds = nan*zeros((len(DDs_use), 2))

    for s, sim in enumerate(sims):
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
        fits_name = '/Users/rsimons/Dropbox/rcs_foggie/satellite_masses/all/%s_all.fits'%(simname)
        #print fits_name
        data = fits.open(fits_name)
        ts_s = data['SAT_%.2i'%sat_n].data['ts']
        ms_s = data['SAT_%.2i'%sat_n].data['ms']
        mg_s = data['SAT_%.2i'%sat_n].data['mg']
        sf_s = data['SAT_%.2i'%sat_n].data['sf']
        dm_s = data['SAT_%.2i'%sat_n].data['dm']

        ax1[0].plot(ts_s, ms_s, label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        ax1[1].plot(ts_s, dm_s, label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        ax1[2].plot(ts_s, mg_s, label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        #axes[1,1].plot(ts, R_90, label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        #ax2.plot(ts_s, sf_s, label = sim, color = clrs[s], linestyle = ls[s], linewidth = 1, alpha = alp[s])

        if s < 4:
            ms_bounds[:,0] = np.nanmin((ms_s, ms_bounds[:,0]), axis = 0)
            dm_bounds[:,0] = np.nanmin((dm_s, dm_bounds[:,0]), axis = 0)
            mg_bounds[:,0] = np.nanmin((mg_s, mg_bounds[:,0]), axis = 0)

            ms_bounds[:,1] = np.nanmax((ms_s, ms_bounds[:,1]), axis = 0)
            dm_bounds[:,1] = np.nanmax((dm_s, dm_bounds[:,1]), axis = 0)
            mg_bounds[:,1] = np.nanmax((mg_s, mg_bounds[:,1]), axis = 0)






    ax1[0].fill_between(ts_s, ms_bounds[:,0], ms_bounds[:,1], color = 'blue', alpha = 0.2)
    ax1[1].fill_between(ts_s, dm_bounds[:,0], dm_bounds[:,1], color = 'blue', alpha = 0.2)
    ax1[2].fill_between(ts_s, mg_bounds[:,0], mg_bounds[:,1], color = 'blue', alpha = 0.2)

    ax1[0].legend(loc = 4)
    ax2.legend(loc = 2)


    fs = 12

    ax1[0].set_ylabel('M$_*$ (M$_{\odot}$)', fontsize = fs)
    ax1[1].set_ylabel('M$_{DM}$ (M$_{\odot}$)', fontsize = fs)
    ax1[2].set_ylabel('M$_{g}$ (M$_{\odot}$)', fontsize = fs)
    ax2.set_ylabel('star formation rate (M$_{\odot}$ yr$^{-1}$)', fontsize = fs)




    for ax in ax1.ravel(): ax.set_xlabel('time (Gyr)', fontsize = fs)
        #ax.set_yscale('log')

    ax2.set_xlabel('time (Gyr)', fontsize = fs)
    fig.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/butterfly/%i_mass.png'%sat_n, dpi = 300)
    #fig2.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/butterfly/%i_SFR.png'%sat_n, dpi = 300)

    plt.close(fig)
    






