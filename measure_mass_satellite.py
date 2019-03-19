import astropy
from astropy.io import fits
import numpy as np
from numpy import *
import math
from joblib import Parallel, delayed
import os, sys, argparse
import yt
import matplotlib.pyplot as plt
import trident

def parse():
    '''
    Parse command line arguments
    ''' 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''\
                                Generate the cameras to use in Sunrise and make projection plots
                                of the data for some of these cameras. Then export the data within
                                the fov to a FITS file in a format that Sunrise understands.
                                ''')
    parser.add_argument('-DD', '--DD', default=None, help='DD to use')

    parser.add_argument('-simname', '--simname', default=None, help='Simulation to be analyzed.')

    parser.add_argument('-haloname', '--haloname', default='halo_008508', help='halo_name')

    args = vars(parser.parse_args())
    return args






if __name__ == '__main__':
    args = parse()
    simname = args['simname']
    DD = int(args['DD'])
    haloname = args['haloname']
    snapname = 'DD%.4i'%DD
    ds = yt.load('/nobackupp2/mpeeples/%s/%s/%s/%s'%(haloname, simname, snapname, snapname))


    def _stars(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 2

    def _youngstars(pfilter, data):
        return data[(pfilter.filtered_type, "age")] < 2.e7

    # these are only the must refine dark matter particles
    def _darkmatter(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 4


    yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
    yt.add_particle_filter("youngstars",function=_youngstars, filtered_type='all',requires=["age"])
    yt.add_particle_filter("darkmatter",function=_darkmatter, filtered_type='all',requires=["particle_type"])

    ds.add_particle_filter('stars')
    ds.add_particle_filter('darkmatter')
    ds.add_particle_filter('youngstars')
    print ('adding trident fields...')
    trident.add_ion_fields(ds, ions=['O VI', 'O VII', 'Mg II', 'Si II', 'C II', 'C III', 'C IV',  'Si III', 'Si IV', 'Ne VIII'])
    
    #Parallel(n_jobs = -1, backend = 'threading')(delayed(measure_mass)(simname = simname, DD = DD, sat_n = sat_n, ds = ds) for sat_n in np.arange(5))
    species_dict = {'H'   : 'H_mass',
                    'H0'   : 'H_p0_mass',
                    'H1'   : 'H_p1_mass',
                    'CII': 'C_p1_mass',
                    'CIII': 'C_p2_mass',
                    'CIV': 'C_p3_mass',
                    'HI': 'H_p0_mass',
                    'MgII': 'Mg_p1_mass',
                    'OVI': 'O_p5_mass',
                    'OVII': 'O_p6_mass',
                    'SiII': "Si_p1_mass",
                    'SiIII': "Si_p2_mass",
                    'SiIV': "Si_p3_mass",
                    'NeVIII': 'Ne_p7_mass',
                    'FeXIV': 'Fe_p13_mass'}

    for sat_n in np.arange(5):
        #cen_file =  np.load('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_sat%.2i_cen.npy'%(simname, DD, sat_n))[()]

        #anchor_xs_box_avg, anchor_ys_box_avg, anchor_zs_box_avg, anchor_vxs_box_avg, anchor_vys_box_avg, anchor_vzs_box_avg, anchor_xs_avg, anchor_ys_avg, anchor_zs_avg, anchor_vxs_avg, anchor_vys_avg, anchor_vzs_avg = cen_file
        #cen = yt.YTArray([anchor_xs_box_avg, anchor_ys_box_avg, anchor_zs_box_avg], 'kpc')


        #cen_np = np.load('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_sat%.2i_1049_cen.npy'%(simname, DD, sat_n))[()]
        cen_fits = fits.open('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_anchorprops.fits'%(simname, DD))
        cenx = cen_fits['SAT_%.2i'%sat_n].data['box_avg'][0]
        ceny = cen_fits['SAT_%.2i'%sat_n].data['box_avg'][1]
        cenz = cen_fits['SAT_%.2i'%sat_n].data['box_avg'][2]
        cen = yt.YTArray([cenx, ceny, cenz], 'kpc')





        central_xyz_fit = np.load('/nobackupp2/rcsimons/foggie_momentum/catalogs/center_%s.npy'%simname)[()]

        xf = central_xyz_fit['x']
        yf = central_xyz_fit['y']
        zf = central_xyz_fit['z']


        r_arr = concatenate((arange(0.25, 2, 0.1), arange(2, 10, 0.25), arange(10, 20, 1.)))

        DM              = []
        gas_tot         = []
        gas_metals      = []
        stars_mass      = []
        stars_youngmass = []
        gas_H           = []
        gas_H0          = []
        gas_H1          = []
        gas_CII         = []
        gas_CIII        = []
        gas_CIV         = []
        gas_OVI         = []
        gas_OVII        = []
        gas_MgII        = []
        gas_SII         = []
        gas_SIII        = []
        gas_SIV         = []
        gas_NeVIII      = []



        '''
        for rr, r in enumerate(r_arr):        
            print (rr, r)
            gc_sphere =  ds.sphere(cen, ds.arr(r,'kpc'))


            DM.append(gc_sphere.quantities.total_quantity([("darkmatter", "particle_mass")]).to('Msun'))
            gas_tot.append(gc_sphere.quantities.total_quantity([("gas", "cell_mass")]).to('Msun'))
            gas_metals.append(gc_sphere.quantities.total_quantity([("gas", "metal_mass")]).to('Msun'))
            stars_mass.append(gc_sphere.quantities.total_quantity([("stars", "particle_mass")]).to('Msun'))
            stars_youngmass.append(gc_sphere.quantities.total_quantity([("youngstars", "particle_mass")]).to('Msun'))
            gas_H.append(gc_sphere.quantities.total_quantity([("gas", species_dict['H'])]).to('Msun'))
            gas_H0.append(gc_sphere.quantities.total_quantity([("gas", species_dict['H0'])]).to('Msun'))
            gas_H1.append(gc_sphere.quantities.total_quantity([("gas", species_dict['H1'])]).to('Msun'))
            gas_CII.append(gc_sphere.quantities.total_quantity([("gas", species_dict['CII'])]).to('Msun'))
            gas_CIII.append(gc_sphere.quantities.total_quantity([("gas", species_dict['CIII'])]).to('Msun'))
            gas_CIV.append(gc_sphere.quantities.total_quantity([("gas", species_dict['CIV'])]).to('Msun'))
            gas_OVI.append(gc_sphere.quantities.total_quantity([("gas", species_dict['OVI'])]).to('Msun'))
            gas_OVII.append(gc_sphere.quantities.total_quantity([("gas", species_dict['OVII'])]).to('Msun'))
            gas_MgII.append(gc_sphere.quantities.total_quantity([("gas", species_dict['MgII'])]).to('Msun'))
            gas_SII.append(gc_sphere.quantities.total_quantity([("gas", species_dict['SiII'])]).to('Msun'))
            gas_SIII.append(gc_sphere.quantities.total_quantity([("gas", species_dict['SiIII'])]).to('Msun'))
            gas_SIV.append(gc_sphere.quantities.total_quantity([("gas", species_dict['SiIV'])]).to('Msun'))
            gas_NeVIII.append(gc_sphere.quantities.total_quantity([("gas", species_dict['NeVIII'])]).to('Msun'))



        #mass = [gas_mass,gas_metal_mass, DM_mass, stars_mass, youngstars_mass]

        #np.save('/nobackupp2/rcsimons/foggie_momentum/satellite_masses/%s_DD%.4i_mass_sat%.2i.npy'%(simname, DD, sat_n), mass)

        master_hdulist = []
        prihdr = fits.Header()
        prihdr['COMMENT'] = "Storing the mass profiles in this FITS file."
        prihdr['simname'] = simname
        prihdr['DDname'] = DD
        prihdr['snapfile'] = '/nobackupp2/mpeeples/%s/%s/%s/%s'%(haloname, simname, snapname, snapname)

        prihdu = fits.PrimaryHDU(header=prihdr)    
        master_hdulist.append(prihdu)

        colhdr = fits.Header()

 

        master_hdulist.append(fits.ImageHDU(data =  array(r_arr         ), header = colhdr, name = 'distance'))
        master_hdulist.append(fits.ImageHDU(data =  array(DM             ), header = colhdr, name = 'dark_matter'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_tot        ), header = colhdr, name = 'gas_tot'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_metals     ), header = colhdr, name = 'gas_metals'))
        master_hdulist.append(fits.ImageHDU(data =  array(stars_mass     ), header = colhdr, name = 'stars_mass'))
        master_hdulist.append(fits.ImageHDU(data =  array(stars_youngmass), header = colhdr, name = 'stars_youngmass'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_H          ), header = colhdr, name = 'gas_H'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_H0         ), header = colhdr, name = 'gas_H0'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_H1         ), header = colhdr, name = 'gas_H1'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_CII        ), header = colhdr, name = 'gas_CII'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_CIII       ), header = colhdr, name = 'gas_CIII'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_CIV        ), header = colhdr, name = 'gas_CIV'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_OVI        ), header = colhdr, name = 'gas_OVI'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_OVII       ), header = colhdr, name = 'gas_OVII'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_MgII       ), header = colhdr, name = 'gas_MgII'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_SII        ), header = colhdr, name = 'gas_SII'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_SIII       ), header = colhdr, name = 'gas_SIII'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_SIV        ), header = colhdr, name = 'gas_SIV'))
        master_hdulist.append(fits.ImageHDU(data =  array(gas_NeVIII     ), header = colhdr, name = 'gas_NeVIII'))
    
        '''


        central_x = xf[0] * DD**4. + xf[1] * DD**3. + xf[2] * DD**2. + xf[3] * DD + xf[4]
        central_y = yf[0] * DD**4. + yf[1] * DD**3. + yf[2] * DD**2. + yf[3] * DD + yf[4]
        central_z = zf[0] * DD**4. + zf[1] * DD**3. + zf[2] * DD**2. + zf[3] * DD + zf[4]


        cen_gal = yt.YTArray([cenx - central_x, ceny - central_y, cenz - central_z], 'kpc')


        R = sqrt(sum(cen_gal**2.))

        sp1 = ds.sphere(cen_gal, ds.arr(R.value + 2.5, 'kpc'))
        sp2 = ds.sphere(cen_gal, ds.arr(R.value - 2.5, 'kpc'))

        shl = ds.intersection([sp1, sp2])
        shl_ad = shl.ds.all_data()  

        shl_dens = ['gas', 'density'].to('g * cm**-3')
        shl_volu = ['index', 'cell_volume'].to('kpc**3.')

        H_dens, edges = np.histogram(np.log10(shl_dens.value), weights = shl_volu.value, bins = np.linspace(-32, -20, 200))

        master_hdulist.append(fits.ImageHDU(data =  array(edges         ), header = colhdr, name = 'gasdens_bins'))
        master_hdulist.append(fits.ImageHDU(data =  array(H_dens             ), header = colhdr, name = 'CGM_gasdens_dist'))

        thdulist = fits.HDUList(master_hdulist)
        fits_name = '/nobackupp2/rcsimons/foggie_momentum/satellite_masses/%s_DD%.4i_mass_sat%.2i.fits'%(simname, DD, sat_n)
        print '\tSaving to ' + fits_name
        thdulist.writeto(fits_name, overwrite = True)


































