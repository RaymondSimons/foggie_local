from astropy.io import fits
plt.ioff()
simnames = ['natural', 'natural_v2', 'natural_v3', 'natural_v4', 'nref11n_nref10f']
anchor_dir = '/Users/rsimons/Dropbox/rcs_foggie/anchor_files'
part_dir = '/Users/rsimons/Dropbox/rcs_foggie/outputs/particles'
DDs = arange(50, 310, 20)

for s, simname in enumerate(simnames):
    for d, DD in enumerate(DDs):
        data = fits.open(part_dir + '/%s_DD%.4i_particles.fits'%(simname, DD))


'''
fig, axes = plt.subplots(1,3, figsize = (10,3))
clrs = ['red', 'blue']
for s, simname in enumerate(simnames):
    for d, DD in enumerate(DDs):
        if '_v' in simname: cenname = simname.replace('natural_v2', 'nref11n_v2_selfshield_z15')
        else: cenname = simname

        anch = fits.open(anchor_dir + '/%s/%s_DD%.4i_anchorprops.fits'%(simname, cenname,  DD))
        ms = anch['SAT_06'].data['anchor_mss']
        gd = anch['sat_06'].data['ids_used_avg']
        gd = gd[gd != 0]
        x = anch['SAT_06'].data['anchor_xs_box'][gd]
        y = anch['SAT_06'].data['anchor_ys_box'][gd]
        z = anch['SAT_06'].data['anchor_zs_box'][gd]


        xmn = anch['SAT_06'].data['box_avg'][0]
        ymn = anch['SAT_06'].data['box_avg'][1]
        zmn = anch['SAT_06'].data['box_avg'][2]

        dx = 50
        for ii, (i, imn) in enumerate(array([(x, xmn), (y, ymn), (z, zmn)])):
            print DD
            axes[ii].hist(i, color = clrs[d], histtype = 'step', bins = arange(imn - dx, imn+ dx, 0.1))
            axes[ii].axvline(imn, linestyle = 'dashed', color = clrs[d], alpha = 0.1)




fig.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/remeasure_centers/test.png', dpi = 300)
'''






















