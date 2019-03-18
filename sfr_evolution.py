import astropy
from astropy.io import fits
import glob
from glob import glob
import scipy
from scipy import interpolate
plt.rcParams['text.usetex'] = True
plt.ioff()
plt.close('all')


selection_DD = 1049


fits_dir = '/Users/rsimons/Dropbox/rcs_foggie/satellite_masses/selected_at_DD%s'%selection_DD
fig_dir = '/Users/rsimons/Dropbox/rcs_foggie/figures/sfr_evolution'


DD_to_t = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/DD_time.npy')[()]
clrs = ['blue', 'red']
prev = 0
#sat_ns_arr = [(0, 0), (4, 2), (3, 4)]
sat_ns_arr = [(4, 2)]

for s, sat_ns in enumerate(sat_ns_arr):
    fig1, ax1 = plt.subplots(1,1, figsize = (6,6))
    fig2, ax2 =  plt.subplots(1,1, figsize = (6,6))

    for sm, simname in enumerate(['nref11n_selfshield_z15', 'nref11n_nref10f_selfshield_z6']):
    #for sm, simname in enumerate(['nref11n_nref10f_selfshield_z6']):
        sat_n = sat_ns[sm]
        DD_arr = arange(700, 1049)
        tt_arr = array([DD_to_t[2][where(DD_to_t[0] == Dnum)] for Dnum in DD_arr])
        sfr_arr = []
        R_90 = []
        for d, DD in enumerate(DD_arr):
            fits_fl = '%s_DD%.4i_mass_sat%.2i_%i.fits'%(simname, DD, sat_n, selection_DD)
            data = fits.open(fits_dir + '/' + fits_fl)
            #print data[6].data[-1]
            frac_mass = data['STARS_MASS'].data[:]/data['STARS_MASS'].data[25]
            f = interpolate.interp1d(frac_mass, data['DISTANCE'].data)
            f2 = interpolate.interp1d(data['DISTANCE'].data, data['STARS_YOUNGMASS'].data/(2.e7))

            R_90.append(f(0.9))
            sfr_arr.append(f2(f(0.9)))



        if simname == 'nref11n_selfshield_z15': lbl = 'natural'
        elif simname == 'nref11n_nref10f_selfshield_z6': lbl = 'forced'
        ax1.plot(tt_arr, array(sfr_arr), '-', color = clrs[sm], label = lbl, linewidth = 3)
        ax2.plot(tt_arr, array(R_90), '-', color = clrs[sm], linewidth = 3)

    ax1.set_ylim(0, 0.8)
    ax1.legend(fontsize = 20)



    #ax.set_yscale('log')
    ax2.set_ylim(0, 3)

    for ax in [ax1, ax2]:
        ax.set_xlabel('time (Gyr)', fontsize = 18)
    ax1.set_ylabel(r'star-formation rate (M$_{\odot}$ yr$^{-1}$)', fontsize = 18)
    ax2.set_ylabel(r'r$_{*,90}$ (kpc)', fontsize = 18)

    figname = fig_dir + '/sat_%.2i.png'%(sat_ns[1])
    fig1.tight_layout()
    fig1.savefig(figname, dpi = 500)
    print 'saved %s'%figname










