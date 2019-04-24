
#sim_name = 'nref11n_nref10f_selfshield_z6'
#sim_name = 'nref11n_selfshield_z15'
import numpy
from numpy import *

#for sim_name in ['natural', 'nref11n_nref10f_selfshield_z6']:
#for sim_name in ['natural', 'nref11n_nref10f_selfshield_z6']:
#for sim_name in ['nref11n_v2_selfshield_z15', 'nref11n_v3_selfshield_z15', 'nref11n_v4_selfshield_z15']:
for sim_name in ['natural']:#, 'nref11n_v2_selfshield_z15', 'nref11n_v3_selfshield_z15', 'nref11n_v4_selfshield_z15']:

    DDmin = 49
    DDmax = 800
    N_split = 9

    sf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/movie/submit_%s_%i_%i_satmass_qsub.sh'%(sim_name, DDmin, DDmax), 'w+')
    for DD in arange(DDmin, DDmax, N_split):
        snap_name = '%.4i_%.4i'%(DD, min(DD + N_split, DDmax))
        sim_snap_name = snap_name + '_' + sim_name +'_satmovie'

        qsub_fname = '%s.qsub'%(sim_snap_name)

        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/movie/%s'%qsub_fname, 'w+')
        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=16:model=has\n')
        qf.write('#PBS -l walltime=1:00:00\n')
        qf.write('#PBS -q normal\n')
        qf.write('#PBS -N %s\n'%sim_snap_name)
        qf.write('#PBS -M rsimons@jhu.edu\n')
        qf.write('#PBS -m abe\n')
        qf.write('#PBS -o ./outfiles/%s_pbs.out\n'%sim_snap_name)
        qf.write('#PBS -e ./outfiles/%s_pbs.err\n'%sim_snap_name)
        qf.write('#PBS -V\n')
        qf.write('#PBS -W group_list=s1938\n\n\n\n')  

        #for DDi in arange(DD, min(DD + N_split, DDmax)):
        out_string =' > ./outfiles/%s_%.4i_%.4i_movie.err > ./outfiles/%s_%.4i_%.4i_movie.out'%(sim_name, DD,min(DD + N_split, DDmax), sim_name, DD,min(DD + N_split, DDmax))
        qf.write('python /u/rcsimons/scripts/foggie_local/make_movie.py -DDmin %i -DDmax %i -simname %s -wd 30 -wdd 6 %s\n'%(DD, min(DD + N_split, DDmax), sim_name, out_string))


        qf.close()


        sf.write('qsub %s\n'%qsub_fname)


    sf.close()