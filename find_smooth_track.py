import matplotlib
matplotlib.use('Agg')
import astropy
from astropy.io import fits
import numpy as np
from numpy import *
import math
from joblib import Parallel, delayed
import os, sys, argparse
import yt
import matplotlib.pyplot as plt


gp_dir = '/Users/rsimons/Dropbox/rcs_foggie/galprops'

simname = 'nref11n_nref10f_selfshield_z6'
for DD in (200, 950):
    a = np.load('%s_DD%.4i_galprops.npy'%(simname, DD))

    