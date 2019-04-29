from astropy.io import fits
simnames = ['natural']
anchor_dir = '/Users/rsimons/Dropbox/rcs_foggie/anchor_files'
DDs = arange(44, 800)

for s, simname in enumerate(simnames):
    for d, DD in enumerate(DDs):
        anch = fits.open(anchor_dir + '/%s/natural_DD%.4i_anchorprops.fits'%(simname, DD))