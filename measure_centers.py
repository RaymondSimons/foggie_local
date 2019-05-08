from astropy.io import fits
plt.ioff()
simnames = ['natural', 'natural_v2', 'natural_v3', 'natural_v4', 'nref11n_nref10f']
#simnames = ['natural']
anchor_dir = '/Users/rsimons/Dropbox/rcs_foggie/anchor_files'
part_dir = '/Users/rsimons/Dropbox/rcs_foggie/outputs/particles'
DDs = arange(50, 190, 20)



fig, axes = plt.subplots(5,3, figsize = (15,25))

for s, simname in enumerate(simnames):
    np.random.seed(1)
    anchors = fits.open('/Users/rsimons/Dropbox/rcs_foggie/figures/select_sats/anchors_%s_DD0150.fits'%simname)

    dis = anchors['SAT02'].data['id']
    dis_rand = np.random.choice(dis[dis > 0], min(100, len(dis)))
    for d, DD in enumerate(DDs):
        data = fits.open(part_dir + '/%s_DD%.4i_particles.fits'%(simname, DD))
        gd = []
        for ii in dis_rand:
            gd_ii = where(data['STARS'].data['id'] == ii)[0]
            if len(gd_ii) > 0:
                gd.append(gd_ii[0])
        gd = array(gd)
        
        print simname, DD, len(gd)

        if len(gd) > 0:
            xs_anchor = data['STARS'].data['x_box'][gd]
            ys_anchor = data['STARS'].data['y_box'][gd]
            zs_anchor = data['STARS'].data['z_box'][gd]
            
            ms_anchor = data['STARS'].data['mass'][gd]
            



            #xbins = arange(13000, 19000, 2)
            axes[s, 0].hist(xs_anchor, weights = ms_anchor, histtype = 'step', color = clrs[s], bins = 200)#xbins)#arange(12450, 12650, 2))
            axes[s, 1].hist(ys_anchor, weights = ms_anchor, histtype = 'step', color = clrs[s], bins = 200)#arange(12450, 12650, 2))
            axes[s, 2].hist(zs_anchor, weights = ms_anchor, histtype = 'step', color = clrs[s], bins = 200)#arange(12500, 12900, 2))


            x = median(xs_anchor)
            y = median(ys_anchor)
            z = median(zs_anchor)

            axes[s, 0].axvline(x, color = clrs[s], linestyle = '--')
            axes[s, 1].axvline(y, color = clrs[s], linestyle = '--')
            axes[s, 2].axvline(z, color = clrs[s], linestyle = '--')



for ax in axes.ravel():
    ax.set_xlim(13825, 19000)

fig.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/select_sats/moving_centers.png', dpi = 300)


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






















