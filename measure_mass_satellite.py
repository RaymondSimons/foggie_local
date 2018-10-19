import matplotlib
matplotlib.use('Agg')
import astropy
from astropy.io import fits
import numpy as np
from numpy import *
import math
from joblib import Parallel, delayed
import os, sys, argparse
import yt
import matplotlib.pyplot as plt

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



def measure_potential(self, r_min = 0.1,  r_step1 = 0.2, r_cen1 = 5, r_step2 = 1,  r_cen2 = 15, r_step3 = 5, r_max = 200.):

    print 'Measuring the potential...'
    center = self.ds.arr([self.cen_x, self.cen_y, self.cen_z], 'kpc')

    rad_steps = concatenate((arange(r_min,  r_cen1, r_step1), 
                             arange(r_cen1, r_cen2, r_step2),
                             arange(r_cen2, r_max,  r_step3)))
    self.mass_profile = zeros((2,len(rad_steps)))

    for i in arange(0,len(rad_steps)):
        print i, rad_steps[i], len(rad_steps)
        try:
            gc_sphere =  self.ds.sphere(center, self.ds.arr(rad_steps[i],'kpc'))
            baryon_mass, particle_mass = gc_sphere.quantities.total_quantity(["cell_mass", "particle_mass"])
            self.mass_profile[0,i] = rad_steps[i]
            self.mass_profile[1,i] = baryon_mass + particle_mass
        except:
            print '\tsomething broken in measure_potential..'
            self.mass_profile[0,i] = 0.
            self.mass_profile[1,i] = 0.


    self.spl = UnivariateSpline(self.mass_profile[0,:], self.mass_profile[1,:])




if __name__ == '__main__':
    args = parse()

    simname = args['simname']
    DD = int(args['DD'])
    haloname = args['haloname']

    snapname = 'DD%.4i'%DD

    ds = yt.load('/nobackupp2/mpeeples/%s/%s/%s/%s'%(haloname, simname, snapname, snapname))

    cen_file =  np.load('/nobackupp2/rcsimons/foggie_momentum/anchor_files/%s_DD%.4i_cen.npy'%(simname, DD))[()]
    anchor_xs_box_avg, anchor_ys_box_avg, anchor_zs_box_avg, anchor_vxs_box_avg, anchor_vys_box_avg, anchor_vzs_box_avg = cen_file

    cen = yt.YTArray([anchor_xs_box_avg, anchor_ys_box_avg, anchor_zs_box_avg], 'kpc')


    gc_sphere =  ds.sphere(center, ds.arr(3,'kpc'))





































