import astropy
from astropy.io import fits
import numpy as np




def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = numpy.average(values, weights=weights)
    # Fast and numerically precise:
    variance = numpy.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))




if __name__ == '__main__':
    min_DD = 610
    max_DD = 611

    momentum_directory = '/nobackupp2/rcsimons/foggie_momentum/momentum_fits'    

    for simname in ['nref11n_selfshield_z15']:
        to_save = []
        anchor_ids = np.load('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_anchors.npy'%anchors)
        for DD in np.arange(min_dd, max_dd):
            fits_name = simname + '_' + 'DD%.4i_momentum.fits'%DD
            ms_s = a['STAR_MASS'].data
            id_s = a['STARS_ID'].data
            x_s, y_s, z_s = a['STARS_GAL_POSITION'].data
            vx_s, vy_s, vz_s = a['STARS_GAL_VELOCITY'].data
            ep_s = a['STARS_EPSILON_FIXED'].data
            x_s_box, y_s_box, z_s_box = a['STARS_BOX_POSITION'].data
            vx_s_box, vy_s_box, vz_s_box = a['STARS_BOX_VELOCITY'].data


            gd_indices = nan * zeros(len(anchor_ids))
            for g in arange(len(anchor_ids)):
                gd_indices[g] = where(id_s = anchor_ids[g])[0]

            anchor_xs    = xs[gd_indices]
            anchor_ys    = ys[gd_indices]
            anchor_zs    = zs[gd_indices]
            anchor_vxs   = vxs[gd_indices]
            anchor_vys   = vys[gd_indices]
            anchor_vzs   = vzs[gd_indices]

            anchor_xs_box    =  xs_box[gd_indices]
            anchor_ys_box    =  ys_box[gd_indices]
            anchor_zs_box    =  zs_box[gd_indices]
            anchor_vxs_box   = vxs_box[gd_indices]
            anchor_vys_box   = vys_box[gd_indices]
            anchor_vzs_box   = vzs_box[gd_indices]


            anchor_xs_box_avg, _  = weighted_avg_and_std(anchor_xs_box, weights = anchor_mss)
            anchor_ys_box_avg, _  = weighted_avg_and_std(anchor_ys_box, weights = anchor_mss)
            anchor_zs_box_avg, _  = weighted_avg_and_std(anchor_zs_box, weights = anchor_mss)
            anchor_vxs_box_avg, _ = weighted_avg_and_std(anchor_vxs_box, weights = anchor_mss)
            anchor_vys_box_avg, _ = weighted_avg_and_std(anchor_vys_box, weights = anchor_mss)
            anchor_vzs_box_avg, _ = weighted_avg_and_std(anchor_vzs_box, weights = anchor_mss)


            to_save.append(['%.4i'%DD, anchor_xs_box_avg, anchor_ys_box_avg, anchor_zs_box_avg, anchor_vxs_box_avg, anchor_vys_box_avg, anchor_vzs_box_avg])


        np.load('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_cens.npy'%to_save)







