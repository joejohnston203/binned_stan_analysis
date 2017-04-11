#!/usr/bin/evn python
#
# Python functions giving the shape of signals and backgrounds
# specific to the Ricochet sensitivity analysis
#
# Each function should take a numpy array of input values, and return
# a numpy array of the output values
#
# ----------------------------------------------------------------------
# Author: J. P. Johnston <jpj13@mit.edu>
#

import numpy as np

def cns_time(time_days):
    if(time_days<365*0.6):
        return 1.0
    else:
        return 0.5
cns_time = np.vectorize(cns_time)

def flat(x):
    return 1.0
flat = np.vectorize(flat)
