import yt
from yt.units import kpc, Mpc
import joblib
from joblib import Parallel, delayed
import seaborn as sns

import os
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from astropy import constants as c
from astropy.cosmology import Planck13 as cosmo
import scipy
from scipy import stats
from yt.fields.api import ValidateParameter

def ram_slice(haloname):
    center_dic =  np.load('/Users/rsimons/Desktop/foggie/ouputs/centers/%s_centers.npy'%haloname, allow_pickle = True)[()]
    to_save = {}

    if haloname == 'halo_008508': 
        simnames =  ['natural', 'nref11c_nref9f', 'nref11n_nref10f'] 
        DDname = 'DD0487'
    else: 
        simnames = ['natural', 'nref11c_nref9f']
        DDname = 'DD0581'

    for simname in simnames[0:1]:
        to_save_sim = {}
        print ('new', haloname, simname)


                
        def _rampressure(field, data):


            cen_bulkv = data.get_field_parameter("bulk_velocity").in_units('km/s')

            velx = data['enzo', 'x-velocity'].to('km/s') - cen_bulkv[0]
            vely = data['enzo', 'y-velocity'].to('km/s') - cen_bulkv[1]
            velz = data['enzo', 'z-velocity'].to('km/s') - cen_bulkv[2]

            vel_mean = np.sqrt(velx**2. + vely**2. + velz**2.)/np.sqrt(3)

            return data['density'] * vel_mean**2.




        yt.add_field(("gas", "ram_pressure"), function=_rampressure, units="kg * km/s**2/kpc**2", validators=[ValidateParameter(['center', 'bulk_velocity'])])
        center = center_dic['center_%s'%simname].to('kpc')

        flname = '/Users/rsimons/Desktop/foggie/sims/%s/%s/%s/%s'%(haloname, simname, DDname, DDname)


        box_size = 100.
        ds = yt.load(flname)
        sp = ds.r[center[0] - yt.YTArray(box_size/2., 'kpc'): center[0] + yt.YTArray(box_size/2., 'kpc'), \
                  center[1] - yt.YTArray(box_size/2., 'kpc'): center[1] + yt.YTArray(box_size/2., 'kpc'), \
                  center[2] - yt.YTArray(box_size/2., 'kpc'): center[2] + yt.YTArray(box_size/2., 'kpc')]

        N = 1
        if not os.path.exists('/Users/rsimons/Dropbox/file_transfer/%s'%haloname):
            os.system('mkdir /Users/rsimons/Dropbox/file_transfer/%s'%haloname)
        if not os.path.exists('/Users/rsimons/Dropbox/file_transfer/%s/%s'%(haloname, simname)):
            os.system('mkdir /Users/rsimons/Dropbox/file_transfer/%s/%s'%(haloname, simname))

        for i in arange(N):
            L = [1*cos(pi*(i)/100.),0, 1*sin(pi*(i)/100.)] # vector normal to cutting plane
            north_vector = [0,1,0]
            num_pix = 512
            #image1 = yt.off_axis_projection(sp, center, L, W, N, ('gas', 'ram_pressure'), north_vector =  north_vector, weight = ('gas', 'cell_mass'))
            #image1 = image1.in_units('kg * km/s**2/kpc**2')
            #fig, axes = plt.subplots(1,1, figsize = (5, 5))
            #im1 = axes.imshow(np.log10(image1), vmin = 18, vmax = 29)
            #fig.savefig('/Users/rsimons/Dropbox/file_transfer/%s/%s/projection_RP_%s_%s_%s.png'%(haloname, simname,haloname, simname, i))



            slc = yt.OffAxisProjectionPlot(ds,  L, 'ram_pressure', center = center, width = (box_size, 'kpc'), weight_field = ('gas', 'cell_mass'), north_vector = north_vector, data_source = sp)
            #slc = yt.OffAxisSlicePlot(ds,  L, 'ram_pressure', center = center, width = (100, 'kpc'))#, north_vector = north_vector)
            pressure_color_map = "Spectral"
            slc.set_zlim(('gas','ram_pressure'), zmin = 1.e19, zmax = 5.e27)
            slc.set_cmap(('gas','ram_pressure'), cmap = pressure_color_map)
            print ('Saving slice...')
            slc.save('/Users/rsimons/Dropbox/file_transfer/%s/%s/slice_RP_%s_%s_%s.png'%(haloname, simname,haloname, simname, i))


            density_proj_min = 5e-2  # msun / pc^2
            density_proj_max = 1e4

            density_color_map = sns.blend_palette(
                ("black", "#4575b4", "#4daf4a", "#ffe34d", "darkorange"), as_cmap=True)

            slc = yt.OffAxisProjectionPlot(ds,  L, 'density', center = center, width = (box_size, 'kpc'), north_vector = north_vector, data_source = sp)
            #slc = yt.OffAxisSlicePlot(ds,  L, 'ram_pressure', center = center, width = (100, 'kpc'))#, north_vector = north_vector)
            #slc.set_zlim(('gas','ram_pressure'), zmin = 1.e19, zmax = 1.e28)
            print ('Saving slice...')
            slc.set_zlim(('gas','density'), zmin = density_proj_min, zmax = density_proj_max)
            slc.set_unit(('gas','density'),'Msun/pc**2')
            slc.set_cmap(('gas','density'),cmap = density_color_map)
            slc.save('/Users/rsimons/Dropbox/file_transfer/%s/%s/slice_density_%s_%s_%s.png'%(haloname, simname,haloname, simname, i))



            '''
            #slc = yt.OffAxisProjectionPlot(ds, L, 'density', center = center, width = (100, 'kpc'))#, north_vector = north_vector)
            slc = yt.OffAxisSlicePlot(ds, L, 'density', center = center, width = (100, 'kpc'))#, north_vector = north_vector)

            print ('Saving slice...')
            slc.save('/Users/rsimons/Dropbox/file_transfer/%s/%s/slice_density_%s_%s_%i.png'%(haloname, simname,haloname, simname, i))
            '''


if __name__ == '__main__':

    halonames = np.array(['halo_002878',
                 'halo_004123',
                 'halo_005016',
                 'halo_005036',
                 'halo_002392',
                 'halo_008508'])
    halonames = np.array([
                 'halo_008508'])






    for haloname in halonames: ram_slice(haloname)



