hdus = []
prim_hdu = fits.PrimaryHDU()
hdus.append(prim_hdu)

for sat_n in arange(2,3):
    anchor_ids = anchor_fits[sat_n]
    gd_indices = []
    for anch_id in anchor_ids: 
        match = where(id_s == anch_id)[0]
        if len(match) > 0: gd_indices.append(int(match))
    gd_indices       = array(gd_indices)

    if len(gd_indices) > 10:
        print 'more than 10 anchor stars found for sat %i in DD%.4i'%(sat_n, DD)
        anchor_mss       = mss[gd_indices]
        anchor_xs_box    =  xs_box[gd_indices]
        anchor_ys_box    =  ys_box[gd_indices]
        anchor_zs_box    =  zs_box[gd_indices]
        anchor_vxs_box   = vxs_box[gd_indices]
        anchor_vys_box   = vys_box[gd_indices]
        anchor_vzs_box   = vzs_box[gd_indices]

        anchor_R = sqrt(anchor_xs_box**2. + anchor_ys_box**2. + anchor_zs_box**2.)

        hist_R, r_edges = np.histogram(anchor_R.value, weights = anchor_mss.value, bins = arange(min(anchor_R.value)-20,max(anchor_R.value)+20, 10))

        Rmid = np.mean([r_edges[argmax(hist_R)], r_edges[argmax(hist_R)+1]])
        good = where(abs(anchor_R.value - Rmid) < 5)[0]


        anchor_xs_box_avg, _  = weighted_avg_and_std(anchor_xs_box,  weights = anchor_mss, good = good)
        anchor_ys_box_avg, _  = weighted_avg_and_std(anchor_ys_box,  weights = anchor_mss, good = good)
        anchor_zs_box_avg, _  = weighted_avg_and_std(anchor_zs_box,  weights = anchor_mss, good = good)
        anchor_vxs_box_avg, _ = weighted_avg_and_std(anchor_vxs_box, weights = anchor_mss, good = good)
        anchor_vys_box_avg, _ = weighted_avg_and_std(anchor_vys_box, weights = anchor_mss, good = good)
        anchor_vzs_box_avg, _ = weighted_avg_and_std(anchor_vzs_box, weights = anchor_mss, good = good)

        box_avg = [anchor_xs_box_avg,
                   anchor_ys_box_avg,
                   anchor_zs_box_avg,
                    anchor_vxs_box_avg,
                    anchor_vys_box_avg,
                    anchor_vzs_box_avg]

        cols1 = fits.ColDefs([fits.Column(name = 'box_avg', array =  box_avg    , format = 'D'),
                              fits.Column(name = 'anchor_mss     ', array =  anchor_mss    , format = 'D'),
                              fits.Column(name = 'anchor_xs_box  ', array =  anchor_xs_box , format = 'D'),
                              fits.Column(name = 'anchor_ys_box  ', array =  anchor_ys_box , format = 'D'),
                              fits.Column(name = 'anchor_zs_box  ', array =  anchor_zs_box , format = 'D'),
                              fits.Column(name = 'anchor_vxs_box ', array =  anchor_vxs_box, format = 'D'),
                              fits.Column(name = 'anchor_vys_box ', array =  anchor_vys_box, format = 'D'),
                              fits.Column(name = 'anchor_vzs_box ', array =  anchor_vzs_box, format = 'D'),
                              fits.Column(name = 'ids_used_avg', array =  good, format = 'I'),
                              ])
    else:
        print 'less than 10 anchor stars found for sat %i in DD%.4i'%(sat_n, DD)
        cols1 = fits.ColDefs([fits.Column(name = 'box_avg     ', array =  None   , format = '0D'),
                              fits.Column(name = 'anchor_mss     ', array =  None, format = '0D'),
                              fits.Column(name = 'anchor_xs_box  ', array =  None, format = '0D'),
                              fits.Column(name = 'anchor_ys_box  ', array =  None, format = '0D'),
                              fits.Column(name = 'anchor_zs_box  ', array =  None, format = '0D'),
                              fits.Column(name = 'anchor_vxs_box ', array =  None, format = '0D'),
                              fits.Column(name = 'anchor_vys_box ', array =  None, format = '0D'),
                              fits.Column(name = 'anchor_vzs_box ', array =  None, format = '0D'),
                              fits.Column(name = 'ids_used_avg', array =  None, format = '0D'),
                              ])

    hdus.append(fits.BinTableHDU.from_columns(cols1, name = 'SAT_%.2i'%sat_n))

hdus_fits = fits.HDUList(hdus)
