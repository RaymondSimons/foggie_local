from astropy.io import fits
import os
import pickle
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
plt.ioff()
plt.close('all')


simnames = ['natural', 
            'natural_v2',
            'natural_v3',
            'natural_v4',
            'nref11n_nref10f',
            'nref11c_nref9f']


simnames = ['natural', 
            'nref11n_nref10f',
            'nref11c_nref9f']


simnames = ['nref11c_nref9f']

#simnames = ['nref11c_nref9f']

DD_arr = arange(50, 1560, 10)
DD_plt = linspace(0, 1800, 1000)
for s, simname in enumerate(simnames):
    print (simname)
    fig, axes = plt.subplots(7,6, figsize = (24, 30))
    coord_datas_all = {}  
    for st in arange(7):
        xs, ys, zs, DDs  = [], [], [], []
        for d, DD in enumerate(DD_arr):
            fname = '/nobackupp2/rcsimons/foggie_momentum/anchor_files/halo_008508/anchor_props/%s/%s_DD%.4i_anchorprops.fits'%(simname, simname, DD)
            if os.path.isfile(fname):
                data = fits.open(fname)
                if st <  6: hd_name = 'SAT_%.2i'%st
                if st == 6: hd_name = 'CENTRAL'
                if len(data[hd_name].data['box_avg']) > 0:
                    DDs.append(DD)
                    xs.append(data[hd_name].data['box_avg'][0])
                    ys.append(data[hd_name].data['box_avg'][1])
                    zs.append(data[hd_name].data['box_avg'][2])



        axes[st, 0].plot(DDs, xs,'ro', markersize = 2, zorder = 1)
        axes[st, 2].plot(DDs, ys,'ro', markersize = 2, zorder = 1)
        axes[st, 4].plot(DDs, zs,'ro', markersize = 2, zorder = 1)



        fx = interp1d(DDs, xs, kind = 7)
        fy = interp1d(DDs, ys, kind = 7)
        fz = interp1d(DDs, zs, kind = 7)

        knd = 3
        knd = 7
        fxe = interp1d(DDs, xs, kind = knd, fill_value = 'extrapolate')
        fye = interp1d(DDs, ys, kind = knd, fill_value = 'extrapolate')
        fze = interp1d(DDs, zs, kind = knd, fill_value = 'extrapolate')



        
        coord_datas = array([(DDs, xs, fx, fxe), (DDs, ys, fy, fye), (DDs, zs, fz, fze)])


        dic_2 = {}
        dic_2['DDs'] = DDs
        dic_2['x'] = xs
        dic_2['fx'] = fx
        dic_2['fxe'] = fxe

        dic_2['y'] = ys
        dic_2['fy'] = fy
        dic_2['fye'] = fye

        dic_2['z'] = zs
        dic_2['fz'] = fz
        dic_2['fze'] = fze


        coord_datas_all[hd_name] = dic_2
        


        for a, coord_data in enumerate(coord_datas):
            DDs_ = coord_data[0]
            xx = coord_data[1]
            f = coord_data[2]
            fe = coord_data[3]
            axes[st,a*2].plot(DDs_, f(DDs_), color = 'blue', linestyle = 'dashed', alpha = 0.3,zorder = 3)
            axes[st,a*2].plot(DD_plt, fe(DD_plt), color = 'black', linestyle = 'dashed', alpha = 0.3,zorder = 3)
            axes[st,a*2+1].plot(DDs_, f(DDs_) - xx, color = 'blue', linestyle = 'dashed', alpha = 0.3,zorder = 3)
            axes[st,a*2+1].plot(DDs_, fe(DDs_) - xx, color = 'black', linestyle = 'dashed', alpha = 0.3,zorder = 3)

        #axes[st,1].plot(DD_plt, fy(DD_plt), color = 'black', linestyle = 'dashed', alpha = 0.3,zorder = 3)
        #axes[st,2].plot(DD_plt, fz(DD_plt), color = 'black', linestyle = 'dashed', alpha = 0.3,zorder = 3)
        if st < 6: ann_str = 'SAT %.2i'%st
        if st == 6: ann_str = 'CENTRAL'
        axes[st,0].annotate(ann_str, (0.1, 0.9), va = 'top', xycoords = 'axes fraction', fontsize = 20)

        for ax in axes[st]: ax.set_xlabel('DD')

        axes[st,0].set_ylabel('x-box (kpc)')
        axes[st,1].set_ylabel('f(DD) - x-box (kpc)')

        axes[st,2].set_ylabel('y-box (kpc)')
        axes[st,3].set_ylabel('f(DD) - y-box (kpc)')

        axes[st,4].set_ylabel('z-box (kpc)')
        axes[st,5].set_ylabel('f(DD) - z-box (kpc)')



    np.save('/nobackupp2/rcsimons/foggie_momentum/catalogs/sat_interpolations/%s_interpolations_DD0150_new.npy'%simname, coord_datas_all, allow_pickle = True)

    #with open('/Users/rsimons/Dropbox/rcs_foggie/outputs/sat_interps/%s_satinterps_DD0150.pkl'%simname, 'wb') as f:
    #    pickle.dump(fx, f)




    fig.savefig('/nobackupp2/rcsimons/foggie_momentum/figures/remeasure_centers/%s_interp_centers_new.png'%simname)





