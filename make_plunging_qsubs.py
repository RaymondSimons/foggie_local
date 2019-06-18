import numpy as np
from numpy import *

sh = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/plunges/run_all.sh', 'w+')
for simname in ['natural', 'natural_v2', 'natural_v3', 'natural_v4', 'nref11n_nref10f', 'nref11c_nref9f']:
    for DD in arange(300, 900, 50):
        snapname = '%s_%i_plunges'%(simname, DD)
        qsub_fname = '%s_%i_plunges.qsub'%(simname, DD)

        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/plunges/%s'%qsub_fname, 'w+')

        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=12:model=san\n')
        qf.write('#PBS -l walltime=0:30:00\n')
        qf.write('#PBS -q normal\n')
        qf.write('#PBS -N %s\n'%snapname)
        qf.write('#PBS -M rsimons@jhu.edu\n')
        qf.write('#PBS -m abe\n')
        qf.write('#PBS -o ./out/%s_pbs.out\n'%snapname)
        qf.write('#PBS -e ./out/%s_pbs.err\n'%snapname)
        qf.write('#PBS -V\n')
        qf.write('#PBS -W group_list=s1698\n\n\n\n')  

        qf.write('python /u/rcsimons/scripts/foggie_local/plunging_orbits.py -simname %s  -DD %i > ./out/%s.err > ./out/%s.out\n'%(simname, DD, snapname, snapname))

        qf.close()

        sh.write('qsub ' + qsub_fname + '\n')