import glob
from glob import glob
from astropy.io import fits
from scipy import interpolate
plt.rcParams['text.usetex'] = True
plt.ioff()
plt.close('all')

DD_to_t = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/DD_time.npy')[()]


DDs_use = arange(100, 1000)

clrs = ['blue', 'navy', 'darkblue', 'royalblue', 'green', 'red']
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
'''
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

'''


tmax = [5.8, 3.1, 2.5, 2.2, 1.6, 1.9]
for sat_n in sats:
    fig = plt. figure(figsize = (15, 10))

    ax_ms = plt.subplot2grid((2, 3), (0,0))
    ax_dm = plt.subplot2grid((2, 3), (0,1))
    ax_mg = plt.subplot2grid((2, 3), (0,2))
    ax_sfr = plt.subplot2grid((2, 3), (1,0), colspan = 3)
    #fig2, ax2 = plt.subplots(1,1, figsize = (12, 6))

    ms_bounds = nan*zeros((len(DDs_use), 2))
    dm_bounds = nan*zeros((len(DDs_use), 2))
    mg_bounds = nan*zeros((len(DDs_use), 2))
    gd = where(ts_s < tmax[sat_n])[0]

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


        ax_ms.plot(ts_s[gd], ms_s[gd], label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        ax_dm.plot(ts_s[gd], dm_s[gd], label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        ax_mg.plot(ts_s[gd], mg_s[gd], label = sim, color = clrs[s], linestyle = ls[s], linewidth = lw[s], alpha = alp[s])
        ax_sfr.plot(ts_s[gd], sf_s[gd], label = sim, color = clrs[s], linestyle = ls[s], linewidth = 1, alpha = alp[s])

        if s < 4:
            ms_bounds[gd,0] = np.nanmin((ms_s[gd], ms_bounds[gd,0]), axis = 0)
            dm_bounds[gd,0] = np.nanmin((dm_s[gd], dm_bounds[gd,0]), axis = 0)
            mg_bounds[gd,0] = np.nanmin((mg_s[gd], mg_bounds[gd,0]), axis = 0)
            ms_bounds[gd,1] = np.nanmax((ms_s[gd], ms_bounds[gd,1]), axis = 0)
            dm_bounds[gd,1] = np.nanmax((dm_s[gd], dm_bounds[gd,1]), axis = 0)
            mg_bounds[gd,1] = np.nanmax((mg_s[gd], mg_bounds[gd,1]), axis = 0)






    ax_ms.fill_between(ts_s[gd], ms_bounds[gd,0], ms_bounds[gd,1], color = 'blue', alpha = 0.2)
    ax_dm.fill_between(ts_s[gd], dm_bounds[gd,0], dm_bounds[gd,1], color = 'blue', alpha = 0.2)
    ax_mg.fill_between(ts_s[gd], mg_bounds[gd,0], mg_bounds[gd,1], color = 'blue', alpha = 0.2)

    ax_ms.legend(loc = 2)
    ax_sfr.legend(loc = 2)


    fs = 15

    ax_ms.set_ylabel('M$_*$ (M$_{\odot}$)', fontsize = fs)
    ax_dm.set_ylabel('M$_{DM}$ (M$_{\odot}$)', fontsize = fs)
    ax_mg.set_ylabel('M$_{g}$ (M$_{\odot}$)', fontsize = fs)
    ax_sfr.set_ylabel('star formation rate (M$_{\odot}$ yr$^{-1}$)', fontsize = fs)




    for ax in [ax_ms, ax_dm, ax_mg, ax_sfr]: 
        ax.set_xlabel('time (Gyr)', fontsize = fs)

    fig.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/butterfly/%i_mass.png'%sat_n, dpi = 300)
    #fig2.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/butterfly/%i_SFR.png'%sat_n, dpi = 300)

    plt.close(fig)
    






