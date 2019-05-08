import astropy
from astropy.io import fits
from scipy.signal import find_peaks
plt.ioff()
simnames = ['natural', 'natural_v2', 'natural_v3', 'natural_v4', 'nref11n_nref10f']
anchor_dir = '/Users/rsimons/Dropbox/rcs_foggie/anchor_files'
part_dir = '/Users/rsimons/Dropbox/rcs_foggie/outputs/particles'
DDmin = 150
DDs = arange(DDmin, DDmin+1, 20)

fig, axes = plt.subplots(5,3, figsize = (15,25))

clrs = ['red', 'green', 'blue', 'purple', 'black']

for s, simname in enumerate(simnames):
    for d, DD in enumerate(DDs):
        data = fits.open(part_dir + '/%s_DD%.4i_particles.fits'%(simname, DD))
        xs, ys, zs = data['STARS'].data['x_box'], data['STARS'].data['y_box'], data['STARS'].data['z_box']
        xd, yd, zd = data['DARK'].data['x_box'], data['DARK'].data['y_box'], data['DARK'].data['z_box']

        ms, md =  data['STARS'].data['mass'], data['DARK'].data['mass']

        ids = data['STARS'].data['id']

        gds = where(xs != 0)
        gdd = where(xd != 0)

        xbins = arange(13825, 14000, 2)

        axes[s, 0].hist(xs[gds], weights = ms[gds], histtype = 'step', color = clrs[s], bins = xbins)#xbins)#arange(12450, 12650, 2))
        axes[s, 1].hist(ys[gds], weights = ms[gds], histtype = 'step', color = clrs[s], bins = 200)#arange(12450, 12650, 2))
        axes[s, 2].hist(zs[gds], weights = ms[gds], histtype = 'step', color = clrs[s], bins = 200)#arange(12500, 12900, 2))


        hst_x, xedges = histogram(xs[gds], weights = ms[gds], bins = xbins)

        peaks, _ = find_peaks(hst_x, height = 0)

        hdus = []
        prim_hdu = fits.PrimaryHDU()
        hdus.append(prim_hdu)


        #by visual inspection, peak = 4 is the central
        for p, peak in enumerate(peaks[4:5]):
            mn_peak = mean([xedges[peak], xedges[peak+1]])
            axes[s,0].axvspan(xmin = mn_peak - 2, xmax = mn_peak + 2, color = 'black', alpha = 0.3)

            gd = where((xs > mn_peak-2) & (xs < mn_peak+2))[0]

            colss = fits.ColDefs([fits.Column(name = 'id'       , array =  data['STARS'].data['id'][gd], format = 'I'),
                                  fits.Column(name = 'mass_DD%.4i'%DDmin    ,  array =  data['STARS'].data['mass'][gd], format = 'D'),
                                  fits.Column(name = 'x_box_DD%.4i'%DDmin ,  array =  data['STARS'].data['x_box'][gd], format = 'D'),
                                  fits.Column(name = 'y_box_DD%.4i'%DDmin ,  array =  data['STARS'].data['y_box'][gd], format = 'D'),
                                  fits.Column(name = 'z_box_DD%.4i'%DDmin ,  array =  data['STARS'].data['z_box'][gd], format = 'D'),
                                  fits.Column(name = 'vx_box_DD%.4i'%DDmin  , array =  data['STARS'].data['vx_box'][gd], format = 'D'),
                                  fits.Column(name = 'vy_box_DD%.4i'%DDmin  , array =  data['STARS'].data['vy_box'][gd], format = 'D'),
                                  fits.Column(name = 'vz_box_DD%.4i'%DDmin  , array =  data['STARS'].data['vz_box'][gd], format = 'D'),
                                  fits.Column(name = 'age_DD%.4i'%DDmin      , array =  data['STARS'].data['age'][gd], format = 'D'),
                                  ])

            hdus.append(fits.BinTableHDU.from_columns(colss, name = 'central'))


        for p, peak in enumerate(peaks):
            mn_peak = mean([xedges[peak], xedges[peak+1]])
            axes[s,0].axvspan(xmin = mn_peak - 2, xmax = mn_peak + 2, color = 'black', alpha = 0.3)

            gd = where((xs > mn_peak-2) & (xs < mn_peak+2))[0]

            colss = fits.ColDefs([fits.Column(name = 'id'       , array =  data['STARS'].data['id'][gd], format = 'I'),
                                  fits.Column(name = 'mass_DD%.4i'%DDmin    ,  array =  data['STARS'].data['mass'][gd], format = 'D'),
                                  fits.Column(name = 'x_box_DD%.4i'%DDmin ,  array =  data['STARS'].data['x_box'][gd], format = 'D'),
                                  fits.Column(name = 'y_box_DD%.4i'%DDmin ,  array =  data['STARS'].data['y_box'][gd], format = 'D'),
                                  fits.Column(name = 'z_box_DD%.4i'%DDmin ,  array =  data['STARS'].data['z_box'][gd], format = 'D'),
                                  fits.Column(name = 'vx_box_DD%.4i'%DDmin  , array =  data['STARS'].data['vx_box'][gd], format = 'D'),
                                  fits.Column(name = 'vy_box_DD%.4i'%DDmin  , array =  data['STARS'].data['vy_box'][gd], format = 'D'),
                                  fits.Column(name = 'vz_box_DD%.4i'%DDmin  , array =  data['STARS'].data['vz_box'][gd], format = 'D'),
                                  fits.Column(name = 'age_DD%.4i'%DDmin      , array =  data['STARS'].data['age'][gd], format = 'D'),
                                  ])

            if p < 4: hdus.append(fits.BinTableHDU.from_columns(colss, name = 'sat%.2i'%p))
            if p > 4: hdus.append(fits.BinTableHDU.from_columns(colss, name = 'sat%.2i'%(p-1)))



        hdus_fits = fits.HDUList(hdus)
        hdus_fits.writeto('/Users/rsimons/Dropbox/rcs_foggie/figures/select_sats/anchors_%s_DD%.4i.fits'%(simname, DDmin), overwrite = True)




for ax in axes.ravel(): 
    ax.set_yscale('log')

for s, simname in enumerate(simnames):
    axes[s,0].annotate(simname, (0.1, 0.9), xycoords = 'axes fraction')


fig.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/select_sats/selecting_anchors.png', dpi = 300)