
#sim_name = 'nref11n_nref10f_selfshield_z6'
#sim_name = 'nref11n_selfshield_z15'
import numpy
from numpy import *

for sim_name in ['nref11n_selfshield_z15', 'nref11n_nref10f_selfshield_z6']:

    ddmin = 660
    ddmax = 1050
    N_split = 5

    sf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/submit_%s_%i_%i_satmass_qsub.sh'%(sim_name, ddmin, ddmax), 'w+')
    for DD in arange(DDmin, DDmax, N_split):
        snap_name = 'DD%.4i_DD%.4i'%(DD, DD + N_split)
        sim_snap_name = snap_name + '_' + sim_name+'_satmass'

        qsub_fname = '%s.qsub'%(sim_snap_name)

        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/%s'%qsub_fname, 'w+')
        
        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=16:model=san\n')
        qf.write('#PBS -l walltime=01:00:00\n')
        qf.write('#PBS -q normal\n')
        qf.write('#PBS -N %s\n'%sim_snap_name)
        qf.write('#PBS -M rsimons@jhu.edu\n')
        qf.write('#PBS -m abe\n')
        qf.write('#PBS -o ./outfiles/%s_pbs.out\n'%sim_snap_name)
        qf.write('#PBS -e ./outfiles/%s_pbs.err\n'%sim_snap_name)
        qf.write('#PBS -V\n')
        qf.write('#PBS -W group_list=s1938\n\n\n\n')  

        for DDi in arange(DD, DD + N_split):
            qf.write('python /u/rcsimons/scripts/foggie_local/measure_mass_satellite.py --DD %i --simname  %s  > ./outfiles/%s.err > ./outfiles/%s.out\n'%(DDi, sim_name, sim_snap_name, sim_snap_name))




        qf.close()





        sf.write('qsub %s\n'%qsub_fname)


    sf.close()