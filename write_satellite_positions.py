import glob
from glob import glob
from astropy.io import fits
import os

anchor_dir = '/Users/rsimons/Dropbox/rcs_foggie/anchor_files'

fls = glob(anchor_dir+'/nref11n_nref10f_selfshield_z6*fits')
fls = sort(fls)

outdir = '/Users/rsimons/for_jesse'
DD_to_t = np.load('/Users/rsimons/Dropbox/rcs_foggie/outputs/DD_time.npy')[()]



fl0  = open(outdir + '/sat00.cat', 'w+')
fl1  = open(outdir + '/sat01.cat', 'w+')
fl2  = open(outdir + '/sat02.cat', 'w+')
fl3  = open(outdir + '/sat03.cat', 'w+')
fl4  = open(outdir + '/sat04.cat', 'w+')
fl5  = open(outdir + '/sat05.cat', 'w+')
fl6  = open(outdir + '/sat06.cat', 'w+')
fl7  = open(outdir + '/sat07.cat', 'w+')
fl8  = open(outdir + '/sat08.cat', 'w+')
fl9  = open(outdir + '/sat09.cat', 'w+')
fl10 = open(outdir + '/sat10.cat', 'w+')

fls_all = array([fl0, 
                 fl1,
                 fl2,
                 fl3,
                 fl4, 
                 fl5, 
                 fl6, 
                 fl7, 
                 fl8, 
                 fl9, 
                 fl10])



for fl in fls_all:
    fl.write('#time redshift x  y z\n')



for DD in arange(700, 1100):
    gd = where(DD_to_t[0] == DD)
    rz = DD_to_t[1][gd][0]
    t = DD_to_t[2][gd][0]

    fl = '/Users/rsimons/Dropbox/rcs_foggie/anchor_files/nref11n_nref10f_selfshield_z6_DD%.4i_anchorprops.fits'%DD
    if os.path.isfile(fl):
        a = fits.open(fl)
        for s, sat_n in enumerate(arange(11)):
            x, y, z = a['SAT_%.2i'%sat_n].data['box_avg'][0:3]
            fls_all[s].write('%.4f \t%.4f \t%.4f \t%.4f \t%.4f \n'%(t, rz, x, y, z))




for fl in fls_all:
    fl.close()



