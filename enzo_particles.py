import sys
import os
import glob
import yt
import numpy as np
from numpy import *
import astropy
from astropy.cosmology import Planck13 as cosmo
import findGalaxyProps as fGP
import os, sys, argparse














if __name__ == "__main__":

    args = parse()

    simname = args['simname']
    snapname = args['snapname']
    haloname = args['haloname']
    run_parallel = args['run_parallel']





