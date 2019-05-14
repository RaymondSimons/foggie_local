import os
from joblib import Parallel, delayed
import numpy as np
from astropy.io import fits
from numpy import *

def retrieve_vdens(simname, DDs):


        if simname == 'natural':
            simfname = 'natural'


        if simname == 'natural_v2':
            simfname = 'nref11n_v2_selfshield_z15'

        
        if simname == 'natural_v3':
            simfname = 'nref11n_v3_selfshield_z15'


        if simname == 'natural_v4':
            simfname = 'nref11n_v4_selfshield_z15'

        
        if simname == 'nref11c_nref9f':
            simfname = 'nref11c_nref9f'

        if simname == 'nref11n_nref10f':
            simfname = 'nref11n_nref10f'


        master_hdulist = []
        prihdr = fits.Header()
        prihdr['COMMENT'] = "Storing the ram pressure percintiles in this FITS file."
        prihdu = fits.PrimaryHDU(header=prihdr)    
        master_hdulist.append(prihdu)
        cols = []
        cols.append(fits.Column(name = 'DDs', array =  np.array(DDs), format = 'D'))
        cols = fits.ColDefs(cols)        
        master_hdulist.append(fits.BinTableHDU.from_columns(cols, name = 'DD'))
        fits_name = '/nobackupp2/rcsimons/foggie_momentum/ram_pressure/percentiles/%s_ram_percentiles.fits'%(simname)


        for sat_n in arange(6):
            dens    = nan*zeros((len(DDs), 5))
            vx      = nan*zeros((len(DDs), 5))
            vy      = nan*zeros((len(DDs), 5))
            vz      = nan*zeros((len(DDs), 5))

            for d, DD in enumerate(DDs):
                data_name = '/nobackupp2/rcsimons/foggie_momentum/ram_pressure/%s/%s_DD%.4i_ram.fits'%(simname, simfname, DD)
                if os.path.isfile(data_name):
                    data = fits.open(data_name)
                    dens[d,:] = np.percentile(data['SAT_%.2i'%sat_n].data['density'], (2, 16, 50, 84, 98)) 
                    vx[d,:]   = np.percentile(data['SAT_%.2i'%sat_n].data['x_velocity'], (2, 16, 50, 84, 98)) 
                    vy[d,:]   = np.percentile(data['SAT_%.2i'%sat_n].data['y_velocity'], (2, 16, 50, 84, 98)) 
                    vz[d,:]   = np.percentile(data['SAT_%.2i'%sat_n].data['z_velocity'], (2, 16, 50, 84, 98)) 
                else:
                    print ('no file found in %s'%data_name)

            cols = []

            for p, perc in enumerate(np.array(['2', '16', '50', '84', '98'])):
                cols.append(fits.Column(name = 'density_%s'%perc,    array =  dens[:,p], format = 'D'))
                cols.append(fits.Column(name = 'x_velocity_%s'%perc, array =  vx[:,p], format = 'D'))
                cols.append(fits.Column(name = 'y_velocity_%s'%perc, array =  vy[:,p], format = 'D'))
                cols.append(fits.Column(name = 'z_velocity_%s'%perc, array =  vz[:,p], format = 'D'))

            cols = fits.ColDefs(cols)
            
            master_hdulist.append(fits.BinTableHDU.from_columns(cols, name = 'SAT_%.2i'%sat_n))
        
        thdulist = fits.HDUList(master_hdulist)
        print ('\tSaving to ' + fits_name)

        thdulist.writeto(fits_name, overwrite = True)







DDs = arange(49, 100)




simnames =  ['natural',
            'natural_v2',
            'natural_v3',
            'natural_v4',
            'nref11c_nref9f',
            'nref11n_nref10f']




Parallel(n_jobs = -1)(delayed(retrieve_vdens)(simname, DDs) for simname in simnames)

























