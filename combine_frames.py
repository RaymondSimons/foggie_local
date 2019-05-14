import os
import numpy as np
from numpy import *


DDs = arange(49, 1000)
simnames =  ['natural',
            'natural_v2',
            'natural_v3',
            'natural_v4',
            'nref11c_nref9f',
            'nref11n_nref10f']


for d, DD in enumerate(DDs):
    for s, simname in enumerate(simnames):
        print (DD, simname)