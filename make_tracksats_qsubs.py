import numpy as np
from numpy import *
split_n = 4

DDmin = 200
DDmax = 1050

for simname in ['nref11n_selfshield_z15', 'nref11n_nref10f_selfshield_z6']:
    sf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/submit_%s_%i_%i_tracksats.sh'%(simname, DDmin, DDmax), 'w+')
    for i in arange(DDmin, DDmax, split_n):
        min_DD = i
        max_DD = i + split_n
        snapname = 'sats_%s_%i_%i'%(simname, min_DD, max_DD)
        qsub_fname = 'tracksats_%s_%i_%i.qsub'%(simname, min_DD, max_DD)
        
        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/%s'%qsub_fname, 'w+')
        
        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=12:model=san\n')
        qf.write('#PBS -l walltime=2:00:00\n')
        qf.write('#PBS -q normal\n')
        qf.write('#PBS -N %s\n'%snapname)
        qf.write('#PBS -M rsimons@jhu.edu\n')
        qf.write('#PBS -m abe\n')
        qf.write('#PBS -o ./outfiles/%s_pbs.out\n'%snapname)
        qf.write('#PBS -e ./outfiles/%s_pbs.err\n'%snapname)
        qf.write('#PBS -V\n')
        qf.write('#PBS -W group_list=s1698\n\n\n\n')  

        qf.write('python /u/rcsimons/scripts/foggie_local/track_satellites.py \
                 -simname %s -DDmin %i -DDmax %i > ./outfiles/%s_track_satellites.err > \
                 ./outfiles/%s_track_satellites.out\n'%(simname, min_DD, max_DD, snapname, snapname))

        qf.close()


        sf.write('qsub %s\n'%qsub_fname)

    sf.close()




