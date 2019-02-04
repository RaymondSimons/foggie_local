import numpy
from numpy import *
import numpy as np
import matplotlib.pyplot as plt


for simname in ['nref11n_nref10f_selfshield_z6', 'nref11n_selfshield_z15']:

    #cen_fits = np.load('/Users/rsimons/Dropbox/rcs_foggie/catalogs/center_%s.npy'%simname)[()]
    cen_fits = np.load('/nobackupp2/rcsimons/foggie_momentum/catalogs/center_%s.npy'%simname)[()]

    xf = cen_fits['x']
    yf = cen_fits['y']
    zf = cen_fits['z']

    #for DDmin in arange(200, 1050, 100):
    DDmin = 200
    DDmax = 1050
    sf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/submit_%s_%i_%i_movie.sh'%(simname, DDmin, DDmax), 'w+')

    for DD_25 in arange(DDmin, DDmax, 25):
        snap_name = 'DD%.4i_DD%.4i'%(DD_25, DD_25 + 25)
        sim_snap_name = snap_name + '_' + simname
        qsub_fname = '%s_movie.qsub'%(sim_snap_name)
        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/%s'%qsub_fname, 'w+')
        
        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=16:model=san\n')
        qf.write('#PBS -l walltime=00:20:00\n')
        qf.write('#PBS -q normal\n')
        qf.write('#PBS -N %s\n'%sim_snap_name)
        qf.write('#PBS -M rsimons@jhu.edu\n')
        qf.write('#PBS -m abe\n')
        qf.write('#PBS -o ./outfiles/%s_pbs.out\n'%sim_snap_name)
        qf.write('#PBS -e ./outfiles/%s_pbs.err\n'%sim_snap_name)
        qf.write('#PBS -V\n')
        qf.write('#PBS -W group_list=s1698\n\n\n\n')      

        for DD in (DD_25, DD_25 + 25):
            cenx = xf[0] * DD**4. + xf[1] * DD**3. + xf[2] * DD**2. + xf[3] * DD + xf[4]
            ceny = yf[0] * DD**4. + yf[1] * DD**3. + yf[2] * DD**2. + yf[3] * DD + yf[4]
            cenz = zf[0] * DD**4. + zf[1] * DD**3. + zf[2] * DD**2. + zf[3] * DD + zf[4]

            out_string =' > ./outfiles/%s_%.4i_movie.err > ./outfiles/%s_%.4i_movie.out'%(simname, DD, simname, DD)
            #qf.write('python /u/rcsimons/scripts/foggie_local/make_movie.py -DD %i -simname %s -cenx %.2f -ceny %.2f -cenz %.2f --lx %.6f --ly %.6f --lz %.6f %s\n'%(DD, simname, cenx, ceny, cenz, Lx, Ly, Lz, out_string))
            qf.write('python /u/rcsimons/scripts/foggie_local/make_movie.py -DD %i -simname %s -cenx %.4f -ceny %.4f -cenz %.4f %s\n'%(DD, simname, cenx, ceny, cenz, out_string))

        qf.close()

        sf.write('qsub %s\n'%qsub_fname)


    sf.close()















