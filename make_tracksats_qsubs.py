import numpy as np
from numpy import *
split_n = 5

DDmin = 44
DDmax = 500

#natural = nref11n_selfshield_z15

for simname in ['natural', 'nref11n_nref10f_selfshield_z6']:
    sf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/tracks/submit_%s_%.4i_%.4i_tracksats.sh'%(simname, DDmin, DDmax), 'w+')
    for i in arange(DDmin, DDmax, split_n):
        min_DD = i
        max_DD = i + split_n
        snapname = 'sats_%s_%.4i_%.4i'%(simname, min_DD, max_DD)
        qsub_fname = 'tracksats_%s_%.4i_%.4i.qsub'%(simname, min_DD, max_DD)
        
        qf = open('/nobackupp2/rcsimons/foggie_momentum/submit_scripts/tracks/%s'%qsub_fname, 'w+')
        
        qf.write('#PBS -S /bin/bash\n')
        qf.write('#PBS -l select=1:ncpus=20:model=ivy\n')
        qf.write('#PBS -l walltime=1:00:00\n')
        qf.write('#PBS -q normal\n')
        qf.write('#PBS -N %s\n'%snapname)
        qf.write('#PBS -M rsimons@jhu.edu\n')
        qf.write('#PBS -m abe\n')
        qf.write('#PBS -o ./outfiles/%s_pbs.out\n'%snapname)
        qf.write('#PBS -e ./outfiles/%s_pbs.err\n'%snapname)
        qf.write('#PBS -V\n')
        qf.write('#PBS -W group_list=s1938\n\n\n\n')  

        qf.write('python /u/rcsimons/scripts/foggie_local/track_satellites.py \
                 -simname %s -DDmin %i -DDmax %i > ./outfiles/%s_track_satellites.err > \
                 ./outfiles/%s_track_satellites.out\n'%(simname, min_DD, max_DD, snapname, snapname))

        qf.close()


        sf.write('qsub %s\n'%qsub_fname)

    sf.close()




