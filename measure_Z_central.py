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
import foggie
from foggie.utils.get_run_loc_etc import get_run_loc_etc
from foggie.utils.get_refine_box import get_refine_box
from foggie.utils.get_halo_center import get_halo_center
from foggie.utils.get_proper_box_size import get_proper_box_size
import os
import argparse
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from foggie.utils.consistency import *
from foggie.utils import yt_fields
from scipy.signal import find_peaks  
import yt
from numpy import *
from photutils.segmentation import detect_sources
from yt.units import kpc
from foggie.utils.foggie_load import *
from astropy.io import ascii

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''identify satellites in FOGGIE simulations''')
    parser.add_argument('-system', '--system', metavar='system', type=str, action='store', \
                        help='Which system are you on? Default is Jase')
    parser.set_defaults(system="pleiades_raymond")
    parser.add_argument('-DD', '--DD', default=None, help='DD to use')

    parser.add_argument('-simname', '--simname', default=None, help='Simulation to be analyzed.')
    parser.add_argument('-simdir', '--simdir', default='/nobackupp2/mpeeples', help='simulation output directory')

    parser.add_argument('-haloname', '--haloname', default='halo_008508', help='halo_name')

    parser.add_argument('--run', metavar='run', type=str, action='store',
                        help='which run? default is natural')
    parser.set_defaults(run="nref11c_nref9f")

    parser.add_argument('--halo', metavar='halo', type=str, action='store',
                        help='which halo? default is 8508 (Tempest)')
    parser.set_defaults(halo="8508")

    parser.add_argument('--pwd', dest='pwd', action='store_true',
                        help='just use the pwd?, default is no')
    parser.set_defaults(pwd=False)

    parser.add_argument('--run_all', dest='run_all', action='store_true',
                        help='just use the pwd?, default is no')
    parser.set_defaults(pwd=False)

    parser.add_argument('--do_sat_proj_plots', dest='do_sat_proj_plots', action='store_true',
                        help='just use the pwd?, default is no')
    parser.set_defaults(pwd=False)

    parser.add_argument('--do_proj_plots', dest='do_proj_plots', action='store_true',
                        help='just use the pwd?, default is no')
    parser.set_defaults(pwd=False)

    parser.add_argument('--do_identify_satellites', dest='do_identify_satellites', action='store_true',
                        help='just use the pwd?, default is no')
    parser.set_defaults(pwd=False)


    parser.add_argument('--do_record_anchor_particles', dest='do_record_anchor_particles', action='store_true',
                        help='just use the pwd?, default is no')
    parser.set_defaults(pwd=False)

    parser.add_argument('--output', metavar='output', type=str, action='store',
                        help='which output? default is RD0020')
    parser.set_defaults(output="RD0020")


    parser.add_argument('--save_dir', metavar='save_dir', type=str, action='store',
                        help='directory to save products')
    parser.set_defaults(save_dir="~")

    parser.add_argument('--use_halo_c_v', dest='use_halo_c_v', action='store_true',
                        help='just use the pwd?, default is no')
    parser.set_defaults(use_halo_c_v=False)

    args = parser.parse_args()
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


if __name__ == '__main__':
    args = parse_args()
    ds, refine_box = load_sim(args)
    gc_sphere =  ds.sphere(ds.halo_center_kpc, ds.arr(20.,'kpc'))
    results = {}
    for star_type in ['stars', 'old_stars', 'young_stars']:
        Z_stars = gc_sphere.quantities.weighted_average_quantity((star_type, 'metallicity_fraction'), weight = (star_type, 'particle_mass')) 
        M_stars = gc_sphere.quantities.total_quantity((star_type, 'particle_mass'))
        results['%s_Z'%star_type] = Z_stars.to('Zsun')
        results['%s_M'%star_type] = M_stars.to('Msun')
        for i in ['x', 'y', 'z']:
            L_stars     = gc_sphere.quantities.total_quantity((star_type, 'particle_angular_momentum_%s'%i))
            results['%s_L_%s'%(gas_type, i)] = L_stars.to('cm**2*g/s')


    for gas_type in ['all_gas', 'cold_gas', 'hot_gas']:
        if gas_type == 'all_gas': low_temp, high_temp = 0, 1.e6
        if gas_type == 'cold_gas': low_temp, high_temp = 0., 1.5e4
        if gas_type == 'hot_gas': low_temp, high_temp = 1.5e4, 1.e6

        temp_cut = gc_sphere.cut_region(["(obj['temperature'] > {}) & (obj['temperature'] < {})".format(low_temp, high_temp)])

        Z_gas = temp_cut.quantities.weighted_average_quantity(('gas', 'metallicity'), weight = ('gas', 'cell_mass')) 
        M_gas = temp_cut.quantities.total_quantity(('gas', 'cell_mass'))
        results['%s_Z'%gas_type] = Z_gas.to('Zsun')
        results['%s_M'%gas_type] = M_gas.to('Msun')
        for i in ['x', 'y', 'z']:
            L_gas     = temp_cut.quantities.total_quantity(('gas', 'angular_momentum_%s'%i))
            results['%s_L_%s'%(gas_type, i)] = L_gas.to('cm**2*g/s')


    fits_name = '/nobackupp2/rcsimons/foggie_momentum/mass_metallicity_momentum/%s_%s_%s_mass.npy'%(args.run, args.halo, args.output)
    np.save(fits_name, results)






























