import numpy as np
from numpy import *
import os, sys, argparse
import glob
import astropy
from astropy.io import fits
import astropy
from astropy.cosmology import Planck15 as cosmo
from joblib import Parallel, delayed
import scipy
from scipy.interpolate import UnivariateSpline



#This file will be used to store the profile of the momentum
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

    parser.add_argument('snap_files', nargs='?', default=None, help='Snapshot files to be analyzed.')

    #parser.add_argument('-s', '--snap_base', default='10MpcBox_csf512_',
    #                    help='Base of the snapshots file names.') 

    #parser.add_argument('-d', '--distance', default=100000, type=float,
    #                    help='Distance between cameras and the center of the galaxy (in [kpc]).')

    #parser.add_argument('--no_export',action='store_true',
    #                    help='Do not export data to fits for Sunrise.') 

    args = vars(parser.parse_args())
    return args





class momentum_obj():
    def __init__(self, simname, aname, snapfile, fits_name):
        self.ds = yt.load(snapfile)
        self.simname = simname
        self.aname = aname
        self.snapfile = snapfile
        self.fits_name = fits_name

    def load(self):
        dd = self.ds.all_data()
        print 'Loading star velocities...'

        try:
            self.stars_vx = dd['io', 'particle_velocity_x'].in_units('km/s')
            assert self.stars_vx.shape > 5
        except AttributeError,AssertionError:
            print "No star particles found, skipping: ", self.ds._file_amr
            return 0


        self.stars_id = dd['io', 'particle_index']
        #self.stars_metallicity1 = dd['io', 'particle_metallicity1']
        #self.stars_metallicity2 = dd['io', 'particle_metallicity2']


        self.stars_vy = dd['io', 'particle_velocity_y'].in_units('km/s')
        self.stars_vz = dd['io', 'particle_velocity_z'].in_units('km/s')

        print 'Loading star positions...'
        self.stars_x = dd['io', 'particle_position_x'].in_units('kpc')
        self.stars_y = dd['io', 'particle_position_y'].in_units('kpc')
        self.stars_z = dd['io', 'particle_position_z'].in_units('kpc')

        print 'Loading star mass...'
        self.star_mass = dd['io', 'particle_mass'].in_units('Msun')

        print 'Loading star age...'
        self.star_creation_time = dd['io', 'creation_time'].in_units('yr')
        self.star_age = self.ds.arr(cosmo.age(self.ds.current_redshift).value, 'Gyr').in_units('yr') - self.star_creation_time
        if False:
            print 'Loading gas velocity...'
            self.gas_vx = dd['gas', 'velocity_x'].in_units('km/s')
            self.gas_vy = dd['gas', 'velocity_y'].in_units('km/s')
            self.gas_vz = dd['gas', 'velocity_z'].in_units('km/s')
            
            print 'Loading gas cell position...'
            self.gas_x = dd['gas', 'x'].in_units('kpc')
            self.gas_y = dd['gas', 'y'].in_units('kpc')
            self.gas_z = dd['gas', 'z'].in_units('kpc')

            print 'Loading gas temperature...'
            self.gas_temp = dd['gas', 'temperature']

            print 'Loading gas cell mass...'
            self.gas_mass = dd['gas', 'cell_mass']


        print 'Finished loading...'
        return 1


    def recenter(self, galprops):
        print 'Recentering...'
        cen_x, cen_y, cen_z = galprops['stars_com'][0]        

        diff = sqrt((amom.stars_x.value - cen_x)**2. + (amom.stars_y.value - cen_y)**2. + (amom.stars_z.value - cen_z)**2.)

        self.id_cen_star      = self.stars_id[argmin(diff)]

        print self.id_cen_star



        '''
        self.cold_cen         = nir_disc_cat[1:4].astype('float')
        self.cen_star_offset  = nir_cat[2:5].astype('float')
        self.cen_star_voffset = nir_cat[5:8].astype('float')
        self.L_disk_s = nir_cat[8:11].astype('float')
        self.L_disk   = nir_disc_cat[7:10].astype('float')
        '''
        '''
        #Determine offset
        self.cen_x  = self.stars_x[self.id_cen_star-1]  - self.ds.arr(self.cen_star_offset[0], 'kpc')
        self.cen_y  = self.stars_y[self.id_cen_star-1]  - self.ds.arr(self.cen_star_offset[1], 'kpc')
        self.cen_z  = self.stars_z[self.id_cen_star-1]  - self.ds.arr(self.cen_star_offset[2], 'kpc')
        self.cen_vx = self.stars_vx[self.id_cen_star-1] - self.ds.arr(self.cen_star_voffset[0], 'km/s')
        self.cen_vy = self.stars_vy[self.id_cen_star-1] - self.ds.arr(self.cen_star_voffset[1], 'km/s')
        self.cen_vz = self.stars_vz[self.id_cen_star-1] - self.ds.arr(self.cen_star_voffset[2], 'km/s')

        #Recenter positions and velocities for stars
        self.stars_x_cen   = self.stars_x  - self.cen_x
        self.stars_y_cen   = self.stars_y  - self.cen_y
        self.stars_z_cen   = self.stars_z  - self.cen_z
        self.stars_pos_cen = array([self.stars_x_cen, self.stars_y_cen, self.stars_z_cen])
        self.stars_pos_mag = sqrt(self.stars_x_cen**2.  + self.stars_y_cen**2.  + self.stars_z_cen**2.)

        self.stars_vx_cen  = self.stars_vx - self.cen_vx
        self.stars_vy_cen  = self.stars_vy - self.cen_vy
        self.stars_vz_cen  = self.stars_vz - self.cen_vz
        self.stars_vel_cen = array([self.stars_vx_cen, self.stars_vy_cen, self.stars_vz_cen])
        self.stars_vel_mag = sqrt(self.stars_vx_cen**2. + self.stars_vy_cen**2. + self.stars_vz_cen**2.)

        #Recenter positions and velocities for gas
        self.gas_x_cen   = self.gas_x  - self.cen_x
        self.gas_y_cen   = self.gas_y  - self.cen_y
        self.gas_z_cen   = self.gas_z  - self.cen_z
        self.gas_pos_cen = array([self.gas_x_cen, self.gas_y_cen, self.gas_z_cen])
        self.gas_pos_mag = sqrt(self.gas_x_cen**2.  + self.gas_y_cen**2.  + self.gas_z_cen**2.)

        self.gas_vx_cen  = self.gas_vx - self.cen_vx
        self.gas_vy_cen  = self.gas_vy - self.cen_vy
        self.gas_vz_cen  = self.gas_vz - self.cen_vz
        self.gas_vel_cen = array([self.gas_vx_cen, self.gas_vy_cen, self.gas_vz_cen])
        self.gas_vel_mag = sqrt(self.gas_vx_cen**2. + self.gas_vy_cen**2. + self.gas_vz_cen**2.)
        '''


    def calc_momentum(self):
        print 'Calculating momentum...'

        #Calculate momentum for stars
        self.stars_jx = self.stars_vz * self.stars_y - self.stars_z * self.stars_vy
        self.stars_jy = self.stars_vx * self.stars_z - self.stars_x * self.stars_vz
        self.stars_jz = self.stars_vy * self.stars_x - self.stars_y * self.stars_vx
        self.stars_j  = array([self.stars_jx, self.stars_jy, self.stars_jz])
        self.stars_j_mag  = sqrt(self.stars_jx**2. + self.stars_jy**2. + self.stars_jz**2.)


        #Calculate momentum for gas
        #self.gas_jx = self.gas_vz * self.gas_y - self.gas_z * self.gas_vy
        #self.gas_jy = self.gas_vx * self.gas_z - self.gas_x * self.gas_vz
        #self.gas_jz = self.gas_vy * self.gas_x - self.gas_y * self.gas_vx
        #self.gas_j  = array([self.gas_jx, self.gas_jy, self.gas_jz])
        #self.gas_j_mag  = sqrt(self.gas_jx**2. + self.gas_jy**2. + self.gas_jz**2.)

    def measure_potential(self, r_min = 0.1,  r_step1 = 0.3, r_cen1 = 3, r_step2 = 2,  r_cen2 = 15, r_step3 = 5, r_max = 45.):

        print 'Measuring the potential...'
        center = self.ds.arr([self.cen_x, self.cen_y, self.cen_z], 'kpc')

        rad_steps = concatenate((arange(r_min,  r_cen1, r_step1), 
                                 arange(r_cen1, r_cen2, r_step2),
                                 arange(r_cen2, r_max,  r_step3)))
        self.mass_profile = zeros((2,len(rad_steps)))

        for i in arange(0,len(rad_steps)):
            print i, rad_steps[i], len(rad_steps)
            gc_sphere =  self.ds.sphere(center, self.ds.arr(rad_steps[i],'kpc'))
            baryon_mass, particle_mass = gc_sphere.quantities.total_quantity(["cell_mass", "particle_mass"])
            self.mass_profile[0,i] = rad_steps[i]
            self.mass_profile[1,i] = baryon_mass + particle_mass
        self.spl = UnivariateSpline(self.mass_profile[0,:], self.mass_profile[1,:])

    def measure_circularity(self):
        print 'Calculating circularity...'

        G = astropy.constants.G.to('kpc^3/Msun*s^2') # in kpc^3/Msun*s^2
        #internal_mass_gas   = self.ds.arr(self.spl(self.gas_pos_mag),'g').in_units('Msun')


        #self.vcirc_gas      = self.ds.arr(sqrt(G*internal_mass_gas/(self.gas_pos_mag)),'kpc/s').in_units('km/s')
        #self.jcirc_gas      = self.vcirc_gas * self.gas_pos_mag

        internal_mass_stars = self.ds.arr(self.spl(self.stars_pos_mag),'g').in_units('Msun')
        self.vcirc_stars    = self.ds.arr(sqrt(G*internal_mass_stars/(self.stars_pos_mag)),'kpc/s').in_units('km/s')
        self.jcirc_stars    = self.vcirc_stars * self.stars_pos_mag
        self.L_mag          = sqrt(self.L_disk[0]**2.+self.L_disk[1]**2.+self.L_disk[2]**2.)
 
        #costheta_gas        = np.dot(self.L_disk, self.gas_j)/(self.gas_j_mag*self.L_mag)
        #self.jz_gas         = costheta_gas*self.gas_j_mag
 
        costheta_stars      = np.dot(self.L_disk, self.stars_j)/(self.stars_j_mag*self.L_mag)
        self.jz_stars       = costheta_stars*self.stars_j_mag
 
        #self.epsilon_gas    = self.jz_gas/self.jcirc_gas
        self.epsilon_stars  = self.jz_stars/self.jcirc_stars



        #costheta_gas   = np.dot(self.L_disk, self.gas_pos)/(self.gas_pos_mag*self.L_mag)
        #self.zz_gas    = self.ds.arr(costheta_gas * self.gas_pos_mag, 'kpc')
        #self.rr_gas    = sqrt(self.gas_pos_mag**2. - self.zz_gas**2.)

        costheta_stars = np.dot(self.L_disk, self.stars_pos)/(self.stars_pos_mag*self.L_mag)
        self.zz_stars  = self.ds.arr(costheta_stars * self.stars_pos_mag, 'kpc')
        self.rr_stars  = sqrt(self.stars_pos_mag**2. - self.zz_stars**2.)

    def gas_momentum_heatmap(self):
        print 'Measuring gas momentum profiles...'

        cold_gas_zz = where((abs(self.rr_gas) < 30) & (self.gas_temp < 1.e4))
       
        eps_min = -2.5
        eps_max = 2.5
        min_z   = -10
        max_z   = 10
        min_r   = 0
        max_r   = 30
        min_rad = 0
        max_rad = 100.
        bins_n  = 200


        '''
        cold_gas_zz = where((abs(self.rr_gas) < max_r) & (self.gas_temp < 1.e4))
        weights = self.gas_mass[cold_gas_zz]

        self.cg_zz_heatmap, self.cg_zz_xedges, self.cg_zz_yedges = np.histogram2d(self.epsilon_gas[cold_gas_zz], self.zz_gas[cold_gas_zz], 
                                                                   bins=[linspace(eps_min,eps_max,bins_n), linspace(min_z,max_z,bins_n)], 
                                                                   weights = weights)


        cold_gas_rr = where((abs(self.zz_gas) < (max_z-min_z)/2.) & (self.gas_temp < 1.e4))
        weights = self.gas_mass[cold_gas_rr]
        print min_r, max_r
        self.cg_rr_heatmap, self.cg_rr_xedges, self.cg_rr_yedges = np.histogram2d(self.epsilon_gas[cold_gas_rr], self.rr_gas[cold_gas_rr], 
                                                                   bins=[linspace(eps_min,eps_max,bins_n), linspace(min_r,max_r,bins_n)], 
                                                                   weights = weights)
        print self.cg_rr_xedges.min()


        cold_gas = where(self.gas_temp < 1.e4)
        weights = self.gas_mass[cold_gas]
        self.cg_rad_heatmap, self.cg_rad_xedges, self.cg_rad_yedges = np.histogram2d(self.epsilon_gas[cold_gas], self.gas_pos_mag[cold_gas], 
                                                                   bins=[linspace(eps_min,eps_max,bins_n), linspace(min_rad,max_rad,bins_n)], 
                                                                   weights = weights)
        '''
    def write_fits(self):
        print '\tGenerating fits for %s...'%self.aname
        master_hdulist = []
        prihdr = fits.Header()
        prihdr['COMMENT'] = "Storing the momentum measurements in this FITS file."
        prihdr['simname'] = self.simname
        prihdr['scale'] = self.aname.strip('a')
        prihdr['snapfile'] = self.snapfile

        prihdu = fits.PrimaryHDU(header=prihdr)    
        master_hdulist.append(prihdu)

        colhdr = fits.Header()

        #master_hdulist.append(fits.ImageHDU(data = nir_mstar_cat                                                            , header = colhdr, name = 'nir_mstar_cat'))
        master_hdulist.append(fits.ImageHDU(data = self.L_disk                                                              , header = colhdr, name = 'nir_net_momentum'))
        #master_hdulist.append(fits.ImageHDU(data = self.L_disk_s                                                            , header = colhdr, name = 'nir_net_momentum_s'))
        master_hdulist.append(fits.ImageHDU(data = self.stars_id                                                            , header = colhdr, name = 'stars_id'))
        #master_hdulist.append(fits.ImageHDU(data = np.stack((self.stars_metallicity1 , self.stars_metallicity2))            , header = colhdr, name = 'stars_metallicity'))
        master_hdulist.append(fits.ImageHDU(data = np.stack((self.stars_x , self.stars_y , self.stars_z))       , header = colhdr, name = 'stars_xyz_position'))
        master_hdulist.append(fits.ImageHDU(data = np.stack((self.stars_vx , self.stars_vy , self.stars_vz))    , header = colhdr, name = 'stars_xyz_velocity'))
        master_hdulist.append(fits.ImageHDU(data = np.stack((self.rr_stars, self.zz_stars))                                 , header = colhdr, name = 'stars_cylindrical_position'))
        master_hdulist.append(fits.ImageHDU(data = np.stack((self.stars_jx, self.stars_jy, self.stars_jz))      , header = colhdr, name = 'stars_xyz_momentum'))
        master_hdulist.append(fits.ImageHDU(data = self.epsilon_stars                                                       , header = colhdr, name = 'stars_epsilon'))
        master_hdulist.append(fits.ImageHDU(data = self.mass_profile                                                        , header = colhdr, name = 'mass_profile'))
        master_hdulist.append(fits.ImageHDU(data = self.star_mass                                                           , header = colhdr, name = 'star_mass'))
        #master_hdulist.append(fits.ImageHDU(data = self.star_creation_time                                                  , header = colhdr, name = 'star_creation_time'))
        master_hdulist.append(fits.ImageHDU(data = self.star_age                                                            , header = colhdr, name = 'star_age'))

        if False:
            master_hdulist.append(fits.ImageHDU(data = np.stack((self.cg_zz_xedges , self.cg_zz_yedges))        , header = colhdr, name = 'gas_zz_epsilon_edges'))
            master_hdulist.append(fits.ImageHDU(data = self.cg_zz_heatmap                                       , header = colhdr, name = 'gas_zz_epsilon'))


            master_hdulist.append(fits.ImageHDU(data = np.stack((self.cg_rr_xedges , self.cg_rr_yedges))        , header = colhdr, name = 'gas_rr_epsilon_edges'))
            master_hdulist.append(fits.ImageHDU(data = self.cg_rr_heatmap                                       , header = colhdr, name = 'gas_rr_epsilon'))


            master_hdulist.append(fits.ImageHDU(data = np.stack((self.cg_rad_xedges , self.cg_rad_yedges))     , header = colhdr, name = 'gas_rad_epsilon_edges'))
            master_hdulist.append(fits.ImageHDU(data = self.cg_rad_heatmap                                     , header = colhdr, name = 'gas_rad_epsilon'))

            master_hdulist.append(fits.ImageHDU(data = np.stack((self.gas_x , self.gas_y , self.gas_z))        , header = colhdr, name = 'gas_xyz_position'))
            master_hdulist.append(fits.ImageHDU(data = np.stack((self.rr_gas, self.zz_gas))                                , header = colhdr, name = 'gas_cylindrical_position'))
            master_hdulist.append(fits.ImageHDU(data = np.stack((self.gas_jx, self.gas_jy, self.gas_jz))       , header = colhdr, name = 'gas_momentum'))
            master_hdulist.append(fits.ImageHDU(data = self.epsilon_gas                                                    , header = colhdr, name = 'gas_epsilon'))
            #master_hdulist.append(fits.ImageHDU(data = self.gas_temp                                                       , header = colhdr, name = 'gas_temperature'))
            #master_hdulist.append(fits.ImageHDU(data = self.gas_mass                                                       , header = colhdr, name = 'gas_mass'))
        

        print '\tSaving to ' + self.fits_name
        thdulist = fits.HDUList(master_hdulist)
        thdulist.writeto(self.fits_name, overwrite = True)

        return master_hdulist


def measure_momentum(snapfile, ds, out_dir):
    mom = None
    print 'Measuring momentum for '+ snapfile
    aname = snapfile.split('/')[-1]
    simname = snapfile.split('/')[-3]
    fits_name = out_dir+'/'+simname+'_'+aname+'_momentum.fits'

    galprops = np.load(simname+'_galprops.npy')
    print 'fits name : ', fits_name

    amom = momentum_obj(simname, aname, snapfile, fits_name)

    check = amom.load()
    if check == 0: return

    '''
    amom.recenter(nir_cat, nir_disc_cat)
    amom.calc_momentum()
    amom.measure_potential()
    amom.measure_circularity()
    amom.gas_momentum_heatmap()
    amom.write_fits(nir_mstar_cat)
    '''
    return amom






if __name__ == "__main__":

    #sargs = parse()
    import yt

    #if args['snap_files'] is not None: snaps = [args['snap_files']]
    #else: 
    #snaps = np.asarray(glob.glob("*.d"))
    

    path_to_snaps = '~/Dropbox/rcs_foggie/data/halo_008508/nref11n_selfshield_z15/RD0018/RD0018'

    #snaps = np.asarray(glob.glob(path_to_snaps))

    snaps = np.asarray(['/Users/rsimons/Dropbox/rcs_foggie/data/halo_008508/nref11n_selfshield_z15/RD0018/RD0018', 
                        '/Users/rsimons/Dropbox/rcs_foggie/data/halo_008508/nref11n_nref10f_selfshield_z6/RD0018/RD0018'])

    print "Generating Sunrise Input for: ", snaps

    #abssnap = os.path.abspath(snaps[0])
    assert os.path.lexists(snaps[0])




    out_dir = '/Users/rsimons/Dropbox/rcs_foggie/outputs'

    assert os.path.lexists(out_dir)



    for i in arange(2):
    
        ds = yt.load(snaps[i])
        ad = ds.all_data()
        #gas_vx = ad['gas', 'velocity_x']
        #amom = measure_momentum(snaps[i], ds, out_dir)
        snapfile = snaps[i]


        print 'Measuring momentum for '+ snapfile
        aname = snapfile.split('/')[-1]
        simname = snapfile.split('/')[-3]
        fits_name = out_dir+'/'+simname+'_'+aname+'_momentum.fits'

        print 'fits name : ', fits_name
        amom = momentum_obj(simname, aname, snapfile, fits_name)

        check = amom.load()
        #if check == 0: return
        galprops = np.load(simname+'_galprops.npy')[()]
        #amom.recenter(galprops)
        cen_x, cen_y, cen_z = galprops['stars_com'][0]        

        diff = sqrt((amom.stars_x.value - cen_x)**2. + (amom.stars_y.value - cen_y)**2. + (amom.stars_z.value - cen_z)**2.)

        amom.id_cen_star      = int(amom.stars_id[argmin(diff)])

        print amom.id_cen_star

        amom.cen_x  = amom.stars_x[argmin(diff)]
        amom.cen_y  = amom.stars_y[argmin(diff)] 
        amom.cen_z  = amom.stars_z[argmin(diff)] 
        amom.cen_vx = amom.stars_vx[argmin(diff)]
        amom.cen_vy = amom.stars_vy[argmin(diff)]
        amom.cen_vz = amom.stars_vz[argmin(diff)]

        amom.stars_x   = amom.stars_x  - amom.cen_x
        amom.stars_y   = amom.stars_y  - amom.cen_y
        amom.stars_z   = amom.stars_z  - amom.cen_z
        amom.stars_pos = array([amom.stars_x, amom.stars_y, amom.stars_z])
        amom.stars_pos_mag = sqrt(amom.stars_x**2.  + amom.stars_y**2.  + amom.stars_z**2.)

        amom.stars_vx  = amom.stars_vx - amom.cen_vx
        amom.stars_vy  = amom.stars_vy - amom.cen_vy
        amom.stars_vz  = amom.stars_vz - amom.cen_vz
        amom.stars_vel = array([amom.stars_vx, amom.stars_vy, amom.stars_vz])
        amom.stars_vel_mag = sqrt(amom.stars_vx**2. + amom.stars_vy**2. + amom.stars_vz**2.)
        if False:
            amom.gas_x   = amom.gas_x  - amom.cen_x
            amom.gas_y   = amom.gas_y  - amom.cen_y
            amom.gas_z   = amom.gas_z  - amom.cen_z
            amom.gas_pos = array([amom.gas_x, amom.gas_y, amom.gas_z])
            amom.gas_pos_mag = sqrt(amom.gas_x**2.  + amom.gas_y**2.  + amom.gas_z**2.)

            amom.gas_vx  = amom.gas_vx - amom.cen_vx
            amom.gas_vy  = amom.gas_vy - amom.cen_vy
            amom.gas_vz  = amom.gas_vz - amom.cen_vz
            amom.gas_vel = array([amom.gas_vx, amom.gas_vy, amom.gas_vz])
            amom.gas_vel_mag = sqrt(amom.gas_vx**2. + amom.gas_vy**2. + amom.gas_vz**2.)



        #amom.L_disk = galprops['gas_L'][0]
        amom.L_disk = [-0.37085436,  0.14802026,  0.91681898]


        amom.calc_momentum()



        amom.measure_potential()


        amom.measure_circularity()
        #amom.gas_momentum_heatmap()


        amom.write_fits()
        '''


        #dirname = os.path.dirname(abssnap)
        #simname = os.path.basename(dirname) #assumes directory name for simulation name



        '''

        '''
        nir_cat_name = simname[0:-2]+'_v2_'+simname[-2:]
        nir_cat = np.loadtxt('/nobackupp2/rcsimons/catalogs/nir_catalogs/GEN3/'+nir_cat_name+'/galaxy_catalogue/Nir_simplified_disc_cat.txt', skiprows = 1, dtype='str')
        nir_disc_cat = np.loadtxt('/nobackupp2/rcsimons/catalogs/nir_catalogs/GEN3/'+nir_cat_name+'/galaxy_catalogue/Nir_disc_cat.txt', skiprows = 1, dtype='str')
        nir_mstar_cat = np.loadtxt('/nobackupp2/rcsimons/catalogs/nir_catalogs/GEN3/'+nir_cat_name+'/galaxy_catalogue/Mstar.txt')


        print "Simulation name:  ", simname

        out_sim_dir = os.path.join('/nobackupp2/rcsimons/momentum_measurements/', simname)
        print os.path.lexists(out_sim_dir)
        if not os.path.lexists(out_sim_dir):
            print 'Creating momentum directory for %s'%out_sim_dir
            os.mkdir(out_sim_dir)                

        new_snapfiles = []
        for sn in snaps:
            aname = sn.split('_')[-1].rstrip('.d')
            snap_dir = os.path.join(simname+'_'+aname+'_sunrise')
            newf = os.path.join(snap_dir,sn)
            new_snapfiles.append(newf)

        new_snapfiles = np.asarray(new_snapfiles)

        #Make Parallel, send 3 at a time to the node (reduce memory overhead)
        Parallel(n_jobs = 2, backend = 'threading')(delayed(measure_momentum)(new_snapfiles[i], out_sim_dir, nir_cat, nir_disc_cat, nir_mstar_cat) for i in arange(len(new_snapfiles)))

        #for i in arange(len(new_snapfiles)):
        #    mom = measure_momentum(new_snapfiles[i], out_sim_dir, nir_cat, nir_disc_cat, nir_mstar_cat)
        '''











