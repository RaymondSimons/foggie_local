#!/u/rcsimons/miniconda3/bin/python3.7
import numpy
from numpy import *




axs = ['x', 'y', 'z']
sats = arange(6)




sf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/movies/submit_all.sh', 'w+')
for sat in sats:
    for a, ax in enumerate(axs):

        snap_name = '%s_%.2i'%(ax, sat)
        sim_snap_name = snap_name + '_movie'

        qsub_fname = '%s.qsub'%(sim_snap_name)

        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/movies/%s'%qsub_fname, 'w+')
        
        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=16:model=ivy\n')
        qf.write('#PBS -l walltime=1:00:00\n')
        qf.write('#PBS -q normal\n')
        qf.write('#PBS -N %s\n'%sim_snap_name)
        qf.write('#PBS -M rsimons@jhu.edu\n')
        qf.write('#PBS -m abe\n')
        qf.write('#PBS -o ./outfiles/%s_pbs.out\n'%sim_snap_name)
        qf.write('#PBS -e ./outfiles/%s_pbs.err\n'%sim_snap_name)
        qf.write('#PBS -V\n')
        qf.write('#PBS -W group_list=s1938\n\n\n\n')  


        qf.write('/u/rcsimons/scripts/foggie_local/combine_frames.py --ax %s --sat  %i > ./outfiles/%s.err > ./outfiles/%s.out\n'%(ax, sat, sim_snap_name, sim_snap_name))


        qf.close()


        sf.write('qsub %s\n'%qsub_fname)


sf.close()