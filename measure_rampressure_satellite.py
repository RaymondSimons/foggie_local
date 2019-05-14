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



def write_ram_fits(ds, cen_name, simname, DD, cen_fits):
        fits_name = '/nobackupp2/rcsimons/foggie_momentum/ram_pressure/%s/%s_DD%.4i_ram.fits'%(cen_name, simname, DD)
        #if os.path.exists(fits_name): return
        if not os.path.isfile(fits_name):
            master_hdulist = []
            prihdr = fits.Header()
            prihdr['COMMENT'] = "Storing the ram pressure values in this FITS file."
            prihdr['simname'] = simname
            prihdr['DDname'] = DD
            prihdr['snapfile'] = '/nobackupp2/mpeeples/%s/%s/%s/%s'%(haloname, simname, DDname, DDname)
            prihdu = fits.PrimaryHDU(header=prihdr)    
            master_hdulist.append(prihdu)
            sat_hdus = []
            for sat_n in np.arange(6):
                hd_name = 'SAT_%.2i'%sat_n


                cenx = cen_fits[hd_name]['fxe'](DD)
                ceny = cen_fits[hd_name]['fye'](DD)
                cenz = cen_fits[hd_name]['fze'](DD)
                cen = yt.YTArray([cenx, ceny, cenz], 'kpc')


                cenx_next = cen_fits[hd_name]['fxe'](DD+10)
                ceny_next = cen_fits[hd_name]['fye'](DD+10)
                cenz_next = cen_fits[hd_name]['fze'](DD+10)


                cen_next = yt.YTArray([cenx_next, ceny_next, cenz_next], 'kpc')

                vec_next =  cen_next - cen
                vec_next = vec_next/np.sqrt(np.sum(vec_next**2.))

                cen_10kpc  = cen + vec_next * yt.YTArray([10], 'kpc')

                cp = ds.cutting(vec_next, cen_10kpc)
                frb = cp.to_frb((15, 'kpc'), 200)

                frb_dens = frb["gas", "density"]
                frb_vx = frb["gas", "velocity_x"].in_units('km/s')
                frb_vy = frb["gas", "velocity_y"].in_units('km/s')
                frb_vz = frb["gas", "velocity_z"].in_units('km/s')
                cols = []
                cols.append(fits.Column(name = 'density', array =  np.array(frb_dens), format = 'D'))
                cols.append(fits.Column(name = 'x_velocity', array =  np.array(frb_vx), format = 'D'))
                cols.append(fits.Column(name = 'y_velocity', array =  np.array(frb_vy), format = 'D'))
                cols.append(fits.Column(name = 'z_velocity', array =  np.array(frb_vz), format = 'D'))
                cols = fits.ColDefs(cols)
                
                sat_hdus.append(fits.BinTableHDU.from_columns(cols, name = hd_name))
                
            master_hdulist.extend(sat_hdus)
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

    ds = yt.load('%s/%s/%s/%s/%s'%(simdir, haloname, simname,  DDname, DDname))

    if simname == 'natural': cen_name = 'natural'
    if 'v2' in simname: cen_name = 'natural_v2'
    if 'v3' in simname: cen_name = 'natural_v3'
    if 'v4' in simname: cen_name = 'natural_v4'
    if simname == 'nref11n_nref10f': cen_name = 'nref11n_nref10f'    
    if simname == 'nref11c_nref9f': cen_name = 'nref11c_nref9f'    

    cen_fits = np.load('/nobackupp2/rcsimons/foggie_momentum/catalogs/sat_interpolations/%s_interpolations_DD0150.npy'%cen_name, allow_pickle=True)[()]

    write_ram_fits(ds, cen_name, simname, DD, cen_fits)























