
#sim_name = 'nref11n_nref10f_selfshield_z6'
#sim_name = 'nref11n_selfshield_z15'
import numpy
from numpy import *

#for sim_name in ['natural', 'nref11n_nref10f_selfshield_z6']:
#for sim_name in ['natural', 'nref11n_nref10f_selfshield_z6']:

simnames = ['natural', 
            'natural_v2',
            'natural_v3',
            'natural_v4',
            'nref11n_nref10f',
            'nref11c_nref9f']

for sim_name in simnames:
    DDmin = 49
    DDmax = 1000
    N_split = 50

    sf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/mass/submit_%s_%i_%i_satmass_qsub.sh'%(sim_name, DDmin, DDmax), 'w+')
    for DD in arange(DDmin, DDmax, N_split):
        snap_name = 'DD%.4i_DD%.4i'%(DD, min(DD + N_split, DDmax))
        sim_snap_name = snap_name + '_' + sim_name+'_satmass'

        qsub_fname = '%s.qsub'%(sim_snap_name)

        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/mass/%s'%qsub_fname, 'w+')
        
        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=16:model=san\n')
        qf.write('#PBS -l walltime=10:00:00\n')
        qf.write('#PBS -q normal\n')
        qf.write('#PBS -N %s\n'%sim_snap_name)
        qf.write('#PBS -M rsimons@jhu.edu\n')
        qf.write('#PBS -m abe\n')
        qf.write('#PBS -o ./outfiles/%s_pbs.out\n'%sim_snap_name)
        qf.write('#PBS -e ./outfiles/%s_pbs.err\n'%sim_snap_name)
        qf.write('#PBS -V\n')
        qf.write('#PBS -W group_list=s1938\n\n\n\n')  

        for DDi in arange(DD, min(DD + N_split, DDmax)):
            qf.write('/u/rcsimons/scripts/foggie_local/measure_mass_satellite.py --DD %i --simname  %s > ./outfiles/%s.err > ./outfiles/%s.out\n'%(DDi, sim_name, sim_snap_name, sim_snap_name))


        qf.close()


        sf.write('qsub %s\n'%qsub_fname)


    sf.close()