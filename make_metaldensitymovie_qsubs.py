import numpy
from numpy import *

for sim_name in ['nref11n_nref10f']:
    DDmin = 44
    DDmax = 1900
    N_split = 20

    sf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/movies/submit_%s_%i_%i_satmass_qsub.sh'%(sim_name, DDmin, DDmax), 'w+')
    for DD in arange(DDmin, DDmax, N_split):
        snap_name = '%.4i_%.4i'%(DD, min(DD + N_split, DDmax))
        sim_snap_name = snap_name + '_' + sim_name +'_satmovie'

        qsub_fname = '%s.qsub'%(sim_snap_name)

        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/movies/%s'%qsub_fname, 'w+')
        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=16:model=san\n')
        qf.write('#PBS -l walltime=5:00:00\n')
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
        qf.write('/u/rcsimons/scripts/foggie_local/make_movie_density_metallicity.py -DDmin %i -DDmax %i -simname %s %s\n'%(DD, min(DD + N_split, DDmax), sim_name, out_string))


        qf.close()


        sf.write('qsub %s\n'%qsub_fname)


    sf.close()

