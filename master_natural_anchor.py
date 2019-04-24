import astropy
from astropy.io import fits
import os


DDs = arange(44, 800)
sat_ns = arange(0, 12)
sims = array(['forced', 'natural', 'natural_v2', 'natural_v3', 'natural_v4'])


arr = nan * zeros((len(sims), len(DDs), len(sat_ns), 3))


anchor_dir = '/Users/rsimons/Dropbox/rcs_foggie/anchor_files'

for d, DD in enumerate(DDs):
    print DD
    for n, sat_n in enumerate(sat_ns):
        for s, sim in enumerate(sims):
            if sim == 'forced':
                dirname = 'nref11n_nref10f_selfshield_z6'
                simname = 'nref11n_nref10f_selfshield_z6'
            if sim == 'natural': 
                dirname = 'natural'
                simname = 'natural'
            if sim == 'natural_v2': 
                dirname = 'natural_v2'
                simname = 'nref11n_v2_selfshield_z15'
            if sim == 'natural_v3': 
                dirname = 'natural_v3'
                simname = 'nref11n_v3_selfshield_z15'
            if sim == 'natural_v4': 
                dirname = 'natural_v4'
                simname = 'nref11n_v4_selfshield_z15'

            anch_file = anchor_dir + '/' + dirname + '/' + simname + '_DD%.4i_anchorprops.fits'%DD
            if os.path.exists(anch_file):
                anch_fits = fits.open(anch_file)
                if ('natutral' in  sim) | (sat_n < 11):
                    box_avg = anch_fits['SAT_%.2i'%sat_n].data['box_avg'][0:3]
                    if len(box_avg) > 0:
                        arr[s, d, n, :] = box_avg
    

 np.save('/Users/rsimons/Dropbox/rcs_foggie/anchor_files/all_anchors.npy', arr)



















 