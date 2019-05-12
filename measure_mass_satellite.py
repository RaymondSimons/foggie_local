#!/u/rcsimons/miniconda3/bin/python3.7
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
import numpy

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
    parser.add_argument('-simdir', '--simdir', default='/nobackupp2/mpeeples', help='simulation output directory')

    parser.add_argument('-haloname', '--haloname', default='halo_008508', help='halo_name')
    parser.add_argument('-parallel', '--parallel',    action = "store_true", default = False, help = 'run in parallel')
    args = vars(parser.parse_args())
    return args



def weighted_quantile(values, quantiles, sample_weight=None, values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of initial array
    :param old_style: if True, will correct output to be consistent with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = numpy.array(values)
    quantiles = numpy.array(quantiles)
    if sample_weight is None:
        sample_weight = numpy.ones(len(values))
    sample_weight = numpy.array(sample_weight)
    assert numpy.all(quantiles >= 0) and numpy.all(quantiles <= 1), 'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = numpy.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = numpy.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= numpy.sum(sample_weight)
    return numpy.interp(quantiles, weighted_quantiles, values)



def write_mass_fits(ds, cen_name, simname, DD, species_dict, species_keys, r_arr, cen_fits):
        fits_name = '/nobackupp2/rcsimons/foggie_momentum/satellite_masses/%s/%s_DD%.4i_mass.fits'%(cen_name, simname, DD)
        if os.path.exists(fits_name): return
        master_hdulist = []
        prihdr = fits.Header()
        prihdr['COMMENT'] = "Storing the mass profiles in this FITS file."
        prihdr['simname'] = simname
        prihdr['DDname'] = DD
        prihdr['snapfile'] = '/nobackupp2/mpeeples/%s/%s/%s/%s'%(haloname, simname, DDname, DDname)
        prihdu = fits.PrimaryHDU(header=prihdr)    
        master_hdulist.append(prihdu)
        for sat_n in np.arange(6):
            cenx = cen_fits['SAT_%.2i'%sat_n]['fxe'](DD)
            ceny = cen_fits['SAT_%.2i'%sat_n]['fye'](DD)
            cenz = cen_fits['SAT_%.2i'%sat_n]['fze'](DD)
            cen = yt.YTArray([cenx, ceny, cenz], 'kpc')
            masses = {}
            for key in species_keys: masses[key] = []
            for rr, r in enumerate(r_arr):        
                print (rr, r)
                print ('Calculating mass inside %i kpc sphere'%r)
                gc_sphere =  ds.sphere(cen, ds.arr(r,'kpc'))
                for key in species_keys: 
                    print (key)
                    masses[key].append(gc_sphere.quantities.total_quantity([species_dict[key]]).to('Msun'))
            cols = []
            cols.append(fits.ImageHDU(name = 'radius', array =  np.array(r_arr), format = 'D'))
            for key in species_keys: 
                cols.append(fits.Column(name = key, array =  np.array(masses[key]), format = 'D'))

            cols = fits.ColDefs(cols)

            master_hdulist.append(fits.BinTableHDU.from_columns(cols, name = 'SAT%.2i'%sat_n))


        thdulist = fits.HDUList(master_hdulist)
        print ('\tSaving to ' + fits_name)

        thdulist.writeto(fits_name, overwrite = True)




if __name__ == '__main__':
    args = parse()
    simname = args['simname'] 
    DD = int(args['DD'])
    haloname = args['haloname']
    DDname = 'DD%.4i'%DD
    simdir = args['simdir']
    in_parallel = args['parallel']

    ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, simname,  DDname, DDname))




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
    species_dict = {'dark_matter'    : ("darkmatter", "particle_mass"),
                    'gas_tot'        : ("gas", "cell_mass"),
                    'gas_metals'     : ("gas", "metal_mass"),
                    'stars_mass'     : ("stars", "particle_mass"),                    
                    'stars_youngmass': ("youngstars", "particle_mass"),
                    'gas_H'      : ("gas", 'H_mass'),
                    'gas_H0'     : ("gas", 'H_p0_mass'),
                    'gas_H1'     : ("gas", 'H_p1_mass'),
                    'gas_CII'    : ("gas", 'C_p1_mass'),
                    'gas_CIII'   : ("gas", 'C_p2_mass'),
                    'gas_CIV'    : ("gas", 'C_p3_mass'),
                    'gas_OVI'    : ("gas", 'O_p5_mass'),
                    'gas_OVII'   : ("gas", 'O_p6_mass'),
                    'gas_MgII'   : ("gas", 'Mg_p1_mass'),
                    'gas_SII'    : ("gas", "Si_p1_mass"),
                    'gas_SIII'   : ("gas", "Si_p2_mass"),
                    'gas_SIV'    : ("gas", "Si_p3_mass"),
                    'gas_NeVIII' : ("gas", 'Ne_p7_mass')}
    # can do the same thing with species_dict.keys(), but it essentially randomizes the order
    species_keys = ['dark_matter',
                    'gas_tot',        
                    'gas_metals',     
                    'stars_mass',     
                    'stars_youngmass',
                    'gas_H',      
                    'gas_H0',     
                    'gas_H1',     
                    'gas_CII',    
                    'gas_CIII',   
                    'gas_CIV',    
                    'gas_OVI',    
                    'gas_OVII',   
                    'gas_MgII',   
                    'gas_SII',    
                    'gas_SIII',   
                    'gas_SIV',    
                    'gas_NeVIII']
    if simname == 'natural': cen_name = 'natural'
    if 'v2' in simname: cen_name = 'natural_v2'
    if 'v3' in simname: cen_name = 'natural_v3'
    if 'v4' in simname: cen_name = 'natural_v4'
    if simname == 'nref11n_nref10f': cen_name = 'nref11n_nref10f'    
    if simname == 'nref11c_nref9f': cen_name = 'nref11c_nref9f'    

    cen_fits = np.load('/nobackupp2/rcsimons/foggie_momentum/catalogs/sat_interpolations/%s_interpolations_DD0150.npy'%cen_name, allow_pickle=True)[()]
    r_arr = np.array([10])
    write_mass_fits(ds, cen_name, simname, DD, species_dict, species_keys, r_arr, cen_fits)























