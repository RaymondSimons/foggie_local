import numpy
from numpy import *
import numpy as np
import matplotlib.pyplot as plt


for simname in ['nref11n_nref10f_selfshield_z6', 'nref11n_selfshield_z15']:
#for simname in ['nref11n_selfshield_z15']:


    #cen_fits = np.load('/Users/rsimons/Dropbox/rcs_foggie/catalogs/center_%s.npy'%simname)[()]
    cen_fits = np.load('/nobackupp2/rcsimons/foggie_momentum/catalogs/center_%s.npy'%simname)[()]

    xf = cen_fits['x']
    yf = cen_fits['y']
    zf = cen_fits['z']

    #for DDmin in arange(200, 1050, 100):
    DDmin = 500
    DDmax = 1050
    N_split = 10
    sf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/submit_%s_%i_%i_movie_sats.sh'%(simname, DDmin, DDmax), 'w+')

    for DD in arange(DDmin, DDmax, N_split):
        snap_name = 'DD%.4i_DD%.4i'%(DD, DD + N_split)
        sim_snap_name = snap_name + '_' + simname
        qsub_fname = '%s_sats_movie.qsub'%(sim_snap_name)
        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/%s'%qsub_fname, 'w+')
        
        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=16:model=has\n')
        qf.write('#PBS -l walltime=0:40:00\n')
        qf.write('#PBS -q normal\n')
        qf.write('#PBS -N %s\n'%sim_snap_name)
        qf.write('#PBS -M rsimons@jhu.edu\n')
        qf.write('#PBS -m abe\n')
        qf.write('#PBS -o ./outfiles/%s_pbs.out\n'%sim_snap_name)
        qf.write('#PBS -e ./outfiles/%s_pbs.err\n'%sim_snap_name)
        qf.write('#PBS -V\n')
        qf.write('#PBS -W group_list=s1698\n\n\n\n')      

        for DDi in arange(DD, DD + N_split):
            #cenx = xf[0] * DDi**4. + xf[1] * DDi**3. + xf[2] * DDi**2. + xf[3] * DDi + xf[4]
            #ceny = yf[0] * DDi**4. + yf[1] * DDi**3. + yf[2] * DDi**2. + yf[3] * DDi + yf[4]
            #cenz = zf[0] * DDi**4. + zf[1] * DDi**3. + zf[2] * DDi**2. + zf[3] * DDi + zf[4]
            sat_n = 0

            out_string =' > ./outfiles/%s_%.4i_movie.err > ./outfiles/%s_%.4i_movie.out'%(simname, DDi, simname, DDi)
            qf.write('python /u/rcsimons/scripts/foggie_local/make_movie.py -DD %i -simname %s -wd 20 -wdd 5 -figname %s_%.4i.png%s\n'%(DDi, simname, simname, DDi, out_string))


        qf.close()

        sf.write('qsub %s\n'%qsub_fname)


    sf.close()















