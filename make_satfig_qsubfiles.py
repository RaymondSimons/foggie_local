
#sim_name = 'nref11n_nref10f_selfshield_z6'
#sim_name = 'nref11n_selfshield_z15'
import numpy
from numpy import *

#for sim_name in ['natural', 'nref11n_nref10f_selfshield_z6']:
#for sim_name in ['natural', 'nref11n_nref10f_selfshield_z6']:
#for sim_name in ['nref11n_v2_selfshield_z15', 'nref11n_v3_selfshield_z15', 'nref11n_v4_selfshield_z15']:
#for sim_name in ['natural', 'nref11n_v2_selfshield_z15', 'nref11n_v3_selfshield_z15', 'nref11n_v4_selfshield_z15']:
for sim_name in ['natural']:#, 'nref11n_v2_selfshield_z15', 'nref11n_v3_selfshield_z15', 'nref11n_v4_selfshield_z15']:

    DDmin = 44
    DDmax = 1570
    N_split = 15

    sf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/movies/submit_%s_%i_%i_satmass_qsub.sh'%(sim_name, DDmin, DDmax), 'w+')
    for DD in arange(DDmin, DDmax, N_split):
        snap_name = '%.4i_%.4i'%(DD, min(DD + N_split, DDmax))
        sim_snap_name = snap_name + '_' + sim_name +'_satmovie'

        qsub_fname = '%s.qsub'%(sim_snap_name)

        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/movies/%s'%qsub_fname, 'w+')
        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=16:model=san\n')
        qf.write('#PBS -l walltime=4:00:00\n')
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
        qf.write('/u/rcsimons/scripts/foggie_local/make_movie_satinterp.py -DDmin %i -DDmax %i -simname %s %s\n'%(DD, min(DD + N_split, DDmax), sim_name, out_string))


        qf.close()


        sf.write('qsub %s\n'%qsub_fname)


    sf.close()




    


DD_temps = [238,
253,
343,
387,
388,
403,
417,
418,
433,
447,
448,
463,
478,
658,
673,
687,
688,
702,
703,
717,
718,
733,
748,
763,
778,
793]
for sim_name in ['nref11n_nref10f', 'nref11c_nref9f']:
    DDmin = 900
    DDmax = 1018
    N_split = 3

    sf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/movies/submit_%s_%i_%i_satmass_qsub.sh'%(sim_name, DDmin, DDmax), 'w+')
    for DD in arange(DDmin, DDmax, N_split):
        snap_name = '%.4i_%.4i'%(DD, min(DD + N_split, DDmax))
        sim_snap_name = snap_name + '_' + sim_name +'_satmovie'

        qsub_fname = '%s.qsub'%(sim_snap_name)

        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/movies/%s'%qsub_fname, 'w+')
        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=16:model=san\n')
        qf.write('#PBS -l walltime=0:30:00\n')
        qf.write('#PBS -q normal\n')
        qf.write('#PBS -N %s\n'%sim_snap_name)
        qf.write('#PBS -M rsimons@jhu.edu\n')
        qf.write('#PBS -m abe\n')
        qf.write('#PBS -o ./outfiles/%s_pbs.out\n'%sim_snap_name)
        qf.write('#PBS -e ./outfiles/%s_pbs.err\n'%sim_snap_name)
        qf.write('#PBS -V\n')
        qf.write('#PBS -W group_list=s1698\n\n\n\n')  

        #for DDi in arange(DD, min(DD + N_split, DDmax)):
        out_string =' > ./outfiles/%s_%.4i_%.4i_movie.err > ./outfiles/%s_%.4i_%.4i_movie.out'%(sim_name, DD,min(DD + N_split, DDmax), sim_name, DD,min(DD + N_split, DDmax))
        qf.write('/u/rcsimons/scripts/foggie_local/make_movie_satinterp.py -DDmin %i -DDmax %i -simname %s %s\n'%(DD, min(DD + N_split, DDmax), sim_name, out_string))


        qf.close()


        sf.write('qsub %s\n'%qsub_fname)


    sf.close()

