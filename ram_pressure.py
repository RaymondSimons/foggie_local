import astropy
from astropy.io import fits
import numpy
from numpy import *
import matplotlib as mpl
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import yt
from IPython.display import Image
from yt.units import G
from astropy.modeling.models import Sersic2D
from astropy.convolution import convolve_fft, Gaussian2DKernel

plt.ioff()
plt.close('all')

snames = ['DD%.4i'%i for i in [906, 907, 956]]

#snames = ['DD%.4i'%i for i in [956, 957]]


#snames = ['DD0906']

fig_dir = '/Users/rsimons/Dropbox/rcs_foggie/figures'
sat_n = 0
for sname in snames:
    W = yt.YTArray([15, 15, 5], 'kpc')
    N = 500
    north_vector = [0,1,0]
    m = fits.open('/Users/rsimons/Dropbox/rcs_foggie/outputs/nref11n_nref10f_%s_momentum.fits'%sname)
    ds = yt.load('/Users/rsimons/Dropbox/rcs_foggie/data/halo_008508/%s/%s'%(sname, sname))
    xs_box, ys_box, zs_box, vxs_box, vys_box, vzs_box, xs, ys, zs, vxs, vys, vzs = \
        np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/nref11n_orig_%s_sat%i.npy'%(sname, sat_n))[()]

    cen   = yt.YTArray([xs_box, ys_box, zs_box], 'kpc')
    L     = yt.YTArray([vxs_box, vys_box, vzs_box], 'km/s')
    L_init = L.copy()
    L_init_mag = sqrt(sum(L_init**2))
    L_mag = sqrt(sum(L**2))
    L_n   = L/L_mag

    x = 16
    fig, ax = plt.subplots(8, x, figsize = ((12/8.)*x,12))
    print 'gas density off-axis projection'
    sat_dens = yt.off_axis_projection(ds, cen, -L_n, W, N, ('gas', 'density'), north_vector =  north_vector)
    print 'gas stellar density off-axis projection'
    sat_mstar_dens = yt.off_axis_projection(ds, cen, -L_n, W, N, ('deposit', 'io_density'), north_vector =  north_vector)

    kern = Gaussian2DKernel(0.5*(N/15.))
    sat_mstar_dens_c = yt.YTArray(convolve_fft(sat_mstar_dens, kern), sat_mstar_dens.units)
    for i in arange(x):
        cenS_i   = cen + yt.YTArray([-x + 2*i], 'kpc') * L_n
        cp       = ds.cutting(L, cenS_i, north_vector)
        frb      = cp.to_frb((15, 'kpc'), N)
        vel_dens = frb["gas", "density"]
        vel_x = L_init[0] - frb["gas", "velocity_x"].in_units('km/s')
        vel_y = L_init[1] - frb["gas", "velocity_y"].in_units('km/s')
        vel_z = L_init[2] - frb["gas", "velocity_z"].in_units('km/s')


        vel_mag  = sqrt(vel_x**2. + vel_y**2. + vel_z**2.)
        vel_dot  = (vel_x * L_n[0] + vel_y * L_n[1] + vel_z * L_n[2])/(vel_mag)
        vel_L    = vel_mag * vel_dot



        if -x + 2*i >= 0:
            ax[0, i].set_title('+%.1f kpc'%(abs(-x + 2*i)), fontsize = 14)
        else:
            ax[0, i].set_title('-%.1f kpc'%(abs(-x + 2*i)), fontsize = 14)

        vmn = median(vel_L.value.ravel()) - 2. * std(vel_L.value.ravel())
        vmx = median(vel_L.value.ravel()) + 2. * std(vel_L.value.ravel())

        ax[0, i].imshow(vel_L.value, vmin = vmn, vmax = vmx, cmap = cm.jet)

        vmn = median(log10(vel_dens.value)) - 2. * std(log10(vel_dens.value))
        vmx = median(log10(vel_dens.value)) + 2. * std(log10(vel_dens.value))
        ax[1, i].imshow(log10(vel_dens.value), vmin = vmn, vmax = vmx, cmap = cm.viridis)

        ram_pres = vel_dens * vel_L**2.
        anch_force = 2 * pi * G * sat_dens * sat_mstar_dens_c

        ax[2, i].imshow(log10(ram_pres.to('dyn/cm**2')), vmin = -20, vmax = -10, cmap = cm.viridis)
        ax[3, i].imshow(np.log10(sat_dens), vmin = -9, vmax = -1)
        ax[4, i].imshow(np.log10(sat_mstar_dens), vmin = -9, vmax = -1)
        ax[5, i].imshow(np.log10(sat_mstar_dens_c), vmin = -9, vmax = -1)
        ax[6,i].imshow(np.log10(anch_force.to('dyn/cm**2')), vmin = -20, vmax = -10, cmap = cm.viridis)
        ax[7,i].imshow(np.log10(ram_pres.to('dyn/cm**2')/anch_force.to('dyn/cm**2')), vmin = -10, vmax = 10, cmap = cm.PiYG)

    fsa = 20
    ax[0,0].annotate('V', (0.7, 0.15), xycoords = 'axes fraction', fontsize= fsa, color = 'white', fontweight = 'bold')
    ax[1,0].annotate(r'$\rho_g$', (0.7, 0.15), xycoords = 'axes fraction', fontsize= fsa, color = 'white', fontweight = 'bold')
    ax[2,0].annotate(r'$P_R$', (0.7, 0.15), xycoords = 'axes fraction', fontsize= fsa, color = 'white', fontweight = 'bold')
    ax[3,0].annotate(r'$\Sigma_g$', (0.7, 0.15), xycoords = 'axes fraction', fontsize= fsa, color = 'white', fontweight = 'bold')
    ax[4,0].annotate(r'$\Sigma_*$', (0.7, 0.15), xycoords = 'axes fraction', fontsize= fsa, color = 'black', fontweight = 'bold')
    ax[5,0].annotate(r'$\Sigma_{*, c}$', (0.7, 0.15), xycoords = 'axes fraction', fontsize= fsa, color = 'white', fontweight = 'bold')
    ax[6,0].annotate(r'$P_A$', (0.7, 0.15), xycoords = 'axes fraction', fontsize= fsa, color = 'white', fontweight = 'bold')
    ax[7,0].annotate(r'$P_R/P_A$', (0.4, 0.15), xycoords = 'axes fraction', fontsize= fsa, color = 'black', fontweight = 'bold')

    for a in ax.ravel():
        a.axis('off')

    fig.subplots_adjust(left = 0.05, right = 0.95, top =0.95, bottom = 0.05, hspace = 0.0, wspace = 0.0)
    fig.savefig(fig_dir + '/density_velocity_map_%s_sat_%i.png'%(sname, sat_n), dpi = 300)











