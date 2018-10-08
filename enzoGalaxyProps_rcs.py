import sys
import os
import glob
import yt
import numpy as np
from numpy import *
import astropy
from astropy.cosmology import Planck13 as cosmo
import findGalaxyProps as fGP
import os, sys, argparse


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

    parser.add_argument('sim_name', nargs='?', default=None, help='Snapshot files to be analyzed.')

    parser.add_argument('snap_name', nargs='?', default=None, help='Snapshot files to be analyzed.')


    args = vars(parser.parse_args())
    return args


if __name__=="__main__":
    #snaps = np.sort(np.asarray(glob.glob("RD????/RD????")))  #ENZO format a list of snapshots in separate directories
    #snaps = np.sort(np.asarray(glob.glob("~/Dropbox/rcs_foggie/data/halo_008508/nref11n_nref10f_selfshield_z6/RD????/RD????")))

    #snaps = np.asarray(['/Users/rsimons/Dropbox/rcs_foggie/data/halo_008508/nref11n_nref10f_selfshield_z6/RD0018/RD0018'])
    form='ENZO'

    args = parse()
    simname = args['sim_name']

    #snaps = np.sort(np.asarray(glob.glob("/nobackupp2/mpeeples/halo_008508/nref11n_selfshield_z15/%s/%s"%(args['snap_name'], args['snap_name']))))
    snaps = np.sort(np.asarray(glob.glob("/nobackupp2/mpeeples/halo_008508/%s/%s/%s"%(args['sim_name'], args['snap_name'], args['snap_name']))))



    assert snaps.shape[0] > 0

    print("Calculating Galaxy Props for "+form+": ", snaps)

    abssnap = os.path.abspath(snaps[0])
    dirname = os.path.dirname(os.path.dirname(abssnap))
    #simname = os.path.basename(dirname) #assumes directory name for simulation name

    print( "Simulation name:  ", simname)
    '''
    particle_headers = []
    particle_data = []
    stars_data = []
    new_snapfiles = []
    for sn in snaps:
        aname=os.path.basename(sn)
        adir=os.path.abspath(os.path.dirname(sn))
        snap_dir = os.path.join(adir,simname+'_'+aname+'_sunrise')
        yt_fig_dir = snap_dir+'/yt_projections'
        print( "Sunrise directory: ", snap_dir)
        if not os.path.lexists(snap_dir):
            print ("Creating Sunrise directory:", snap_dir)
            os.mkdir(snap_dir)        
        if not os.path.lexists(yt_fig_dir):
            print ("Creating YT figure directory:", yt_fig_dir)
            os.mkdir(yt_fig_dir)        


        new_snapfiles.append(os.path.abspath(sn))
    new_snapfiles = np.asarray(new_snapfiles)
    '''

    new_snapfiles = np.asarray(snaps)


    galaxy_props = {}
    fields = ['scale', 'stars_total_mass', 'stars_com', 'stars_maxdens', 'stars_maxndens', 'stars_hist_center',
              'stars_rhalf', 'stars_mass_profile', 'stars_L','gas_total_mass', 'gas_maxdens', 'gas_L', 'rvir', 
              'Mvir_dm', 'stars_center','snap_files']
    for field in fields: 
        if field in ['scale', 'stars_total_mass', 'stars_rhalf', 'gas_total_mass' ]:
            galaxy_props[field] = np.array([])                
        else:
            galaxy_props[field] = []



    def _stars(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 2

    #this gets dark matter particles in zoom region only
    def _darkmatter(pfilter, data):
        return data[(pfilter.filtered_type, "particle_type")] == 4

    yt.add_particle_filter("stars",function=_stars, filtered_type='all',requires=["particle_type"])
    yt.add_particle_filter("darkmatter",function=_darkmatter, filtered_type='all',requires=["particle_type"])

    ts = yt.DatasetSeries(new_snapfiles)


    for ds,snap_dir in zip(reversed(ts),np.flipud(new_snapfiles)):
        print( "Getting galaxy props: ",  snap_dir)

        ds.add_particle_filter('stars')
        ds.add_particle_filter('darkmatter')

        dd = ds.all_data()
        ds.domain_right_edge = ds.arr(ds.domain_right_edge,'code_length')
        ds.domain_left_edge  = ds.arr(ds.domain_left_edge,'code_length')

        try:
            print 'Loading data...'
            stars_pos_x = dd['stars', 'particle_position_x'].in_units('kpc')
            print 'Loaded.'
            assert stars_pos_x.shape[0] > 5
        except AttributeError:
            print("No star particles found, skipping: ", snap_dir)
            continue


        scale = round(1.0/(ds.current_redshift+1.0),3)
        galaxy_props['scale'] = np.append(galaxy_props['scale'], scale)
    
        galaxy_props['snap_files'] = np.append(galaxy_props['snap_files'],snap_dir)


        print( 'Determining center...')
        max_ndens_arr = fGP.find_center(dd, ds, cen_pos = ds.domain_center.in_units('kpc')[0].value[()], units = 'kpc')
        print( '\tCenter = ', max_ndens_arr)
        sys.stdout.flush()

        #Generate Sphere Selection
        print( 'Determining virial radius...')
        rvir = fGP.find_rvirial(dd, ds, max_ndens_arr)
        print( '\tRvir = ', rvir)
        sys.stdout.flush()

        hc_sphere = ds.sphere(max_ndens_arr, rvir)

 
        galaxy_props['stars_maxndens'].append(max_ndens_arr.value)
        galaxy_props['rvir'] = np.append(galaxy_props['rvir'], rvir.value[()])
        galaxy_props['Mvir_dm'] = np.append(galaxy_props['Mvir_dm'], hc_sphere[('darkmatter', 'particle_mass')].in_units('Msun').sum().value[()])

        
        #Find Galaxy Properties
        galaxy_props = fGP.find_galaxyprops(galaxy_props, ds, hc_sphere, max_ndens_arr)


        del (hc_sphere)
        sys.stdout.flush()




    # Save galaxy props file
    galprops_outdir = '/nobackupp2/rcsimons/foggie_momentum/galprops'
    galaxy_props_file = galprops_outdir + '/' + simname + '_' + args['snap_name'] + '_galprops.npy'


    print( '\nSuccessfully computed galaxy properties')
    print( 'Saving galaxy properties to ', galaxy_props_file)

    np.save(galaxy_props_file, galaxy_props)  



















