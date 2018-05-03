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
import ROOT
import random
import math

import scipy.integrate as spint

# Functions for time dependence
# ----------------------------------------------------------------------
# Uniform time dependence
def flat(x,params=[]):
    return 1.0
flat = np.vectorize(flat)
# ----------------------------------------------------------------------
# falling_exp:
# Falling exponential with decay constant alpha, so y=e^(-alpha*x).
# The default decay constant is 0, which gives a flat time dependence
def falling_exp_0(x,alpha):
    return np.exp(-alpha*x)
falling_exp_0 = np.vectorize(falling_exp_0)
def falling_exp(x,params):
    if(len(params)>0 and params[0]!=""):
        alpha = params[0]
    else:
        alpha = 0.0
    return falling_exp_0(x,alpha)
# ----------------------------------------------------------------------
def gaus_norm(x,mu,sigma):
    return 1.0/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/2.0/sigma**2)
def falling_exp_convolv_temp(x,mu,sigma,lb,ub,params):
    return spint.quad(lambda t:falling_exp(x-t,params)*gaus_norm(t,mu,sigma),lb,ub)[0]
def falling_exp_convolved(x,params):
    return falling_exp_convolv_temp(x,0.0,15,-100,100,params)
falling_exp_convolved = np.vectorize(falling_exp_convolved)
# ----------------------------------------------------------------------
# rate_cns_shape:
# Rate shape at Double Chooz, defined s.t. the value at full power is 1. 
r0 = 400.0
r1a = 355.4
r1b = 468.8
ran_cns_shape = ROOT.TRandom3()
ran_cns_shape.SetSeed(0)
def rate_cns_shape(mode,fluctuation_frac=0.05):
    if(mode==0):
        rate = 0.0
    elif(mode==1):
        rate = 0.95
    elif(mode==2):
        rate = 0.5*(r0/r1a)**2
    elif(mode==3):
        rate = 0.5*(r0/r1b)**2
    else:
        print("Error- invalid mode for rate_cns_shape in data_generator_ricochet.py")
    if(fluctuation_frac>0.0):
        return max(0.0,ran_cns_shape.Gaus(rate,fluctuation_frac*rate))
    else:
        return rate
rate_cns_shape = np.vectorize(rate_cns_shape)
# ----------------------------------------------------------------------
# cns_time:
# Gives the rate as a function of time at Double Chooz, normalized
# such that the maximum value is 1.0. The rate in each bin is randomly
# chosen to be for both reactors on, one reactor on, or both off. Bins
# can be correlated so the reactors tend to stay in the same
# configuration for a given number of days
# parameters:
#   params[0]: full_power_frac = frac of the time with both reactors on
#   params[1]: off_frac = frac of the time with both reactors off
#     For the remaining time, one of the reactors is on
#   params[2]: force_fractions = If true, then the time spent in each
#              state is forced to be the given fraction of the total
#              time (Within rounding error, since nbins must be integer).
#              (Default = True)
#   params[3]: correlate_bins = If true, then the reactor will tend to
#              stay in the same state
#              (Default = True)
#   params[4]: correlate_days = Number of days the reactor will tend to
#              stay in the same state if correlate_bins==True
#              (Default = 40)
#   params[5]: fluctuation_frac = fraction that the power is allowed to
#              gaussian fluctuate by in each bin.
def cns_time_0(time_days, full_power_frac=0.6,off_frac=0.0,
               force_fractions=True,correlate_bins=True,correlate_days=40,
               fluctuation_frac=0.05):
    ran = ROOT.TRandom3()
    ran.SetSeed(0)
    nbins = len(time_days)
    full_bins = int(full_power_frac*nbins)
    off_bins = int(off_frac*nbins)
    half_1_bins = (nbins-full_bins-off_bins)/2
    half_2_bins = nbins-full_bins-off_bins-half_1_bins
    if(correlate_bins):
        #duration = time_days[-1]-time_days[0]
        #bin_width = duration/float(nbins)
        # Assume equally spaced bins
        bin_width = time_days[1]-time_days[0]
        # An actual number of bins will be poisson distributed around
        # expected_bins
        expected_bins = float(correlate_days)/bin_width
        if(expected_bins==0.0):
            correlate_bins=False
    output_shape = np.zeros(nbins)
    i=0
    while i<nbins:
        if(correlate_bins):
            bins = max(1,ran.Poisson(expected_bins))
        else:
            bins=1
        x = ran.Uniform()
        if(x<full_power_frac):
            # Both reactors on
            mode = 1
            full_bins-=bins
            if(force_fractions and full_bins<0):
                bins+=full_bins
                full_bins=0
        elif(x<full_power_frac+off_frac):
            # Both reactors off
            mode = 0
            off_bins-=bins
            if(force_fractions and off_bins<0):
                bins+=off_bins
                off_bins=0
        elif(x<(1.0+full_power_frac+off_frac)/2.0):
            # Only reactor 1 for half the remaining time
            mode = 2
            half_1_bins-=bins
            if(force_fractions and half_1_bins<0):
                bins+=half_1_bins
                half_1_bins=0
        else:
            # Only reactor 2 for the remaining time
            mode = 3
            half_2_bins-=bins
            if(force_fractions and half_2_bins<0):
                bins+=half_2_bins
                half_2_bins=0
        while(bins>0 and i<nbins):
            output_shape[i] = rate_cns_shape(mode,fluctuation_frac)
            i+=1
            bins-=1
    return output_shape
def cns_time(time_days,params):
    if(len(params)>0 and params[0]!=""):
        full_power_frac = params[0]
    else:
        full_power_frac = 0.6
    if(len(params)>1 and params[1]!=""):
        off_frac = params[1]
    else:
        off_frac = 0.0
    if(len(params)>2 and params[2]!=""):
        force_fractions = params[2]
    else:
        force_fractions = True
    if(len(params)>3 and params[3]!=""):
        correlate_bins=params[3]
    else:
        correlate_bins=True
    if(len(params)>4 and params[4]!=""):
        correlate_days=params[4]
    else:
        correlate_days=40
    if(len(params)>5 and params[5]!=""):
        fluctuation_frac=params[5]
    else:
        fluctuation_frac=0.05
    return cns_time_0(time_days,full_power_frac,off_frac,force_fractions,\
                      correlate_bins,correlate_days,fluctuation_frac)
# ----------------------------------------------------------------------
# cns_time_onillon:
# Approximate time dependence of neutrino signal in Onillon thesis,
# used for making a nice plot. Gives 60\% both on, 0% both off
def time_profile_onillon(x):
    if(0.0<=x<10.0/365.0):
        return 1
    elif(x<20.0/365.0):
        return 2
    elif(x<90.0/365.0):
        return 1
    elif(x<100.0/365.0):
        return 3
    elif(x<120.0/365.0):
        return 1
    elif(x<125.0/365.0):
        return 2
    elif(x<170.0/365.0):
        return 1
    elif(x<228.0/365.0):
        return 2
    elif(x<302.0/365.0):
        return 1
    elif(x<=1.0):
        return 3
    else:
        print("time frac > 1. Returning full power")
        return 1
def cns_time_onillon(time_days,params=[]):
    while time_days>365.0:
        time_days-=365.0
    return rate_cns_shape(time_profile_onillon(time_days/365.0))
cns_time_onillon = np.vectorize(cns_time_onillon)
# ----------------------------------------------------------------------
# cns_time_simple:
# Function with simple time dependence where every third of the year
# first has the given fraction of full power bins, then the given
# fraction of bins with both off, then the rest of the third has
# one reactor on
def cns_time_simple_0(time_days,full_power_frac,off_frac,fluctuation_frac):
    time_days*=3
    while time_days>365.0:
        time_days-=365.0
    x = time_days/365.0
    if(x<full_power_frac):
        # Both reactors on
        mode = 1
    elif(x<full_power_frac+off_frac):
        # Both reactors off
        mode = 0
    elif(x<(1.0+full_power_frac+off_frac)/2.0):
        # Only reactor 1 for half the remaining time
        mode = 2
    else:
        # Only reactor 2 for the remaining time
        mode = 3
    return rate_cns_shape(mode,fluctuation_frac)
cns_time_simple_0 = np.vectorize(cns_time_simple_0)
# params is a list containing up to two elements (enter "" for an
# element to use the Default value):
#   params[0] = fraction of the time at full power (Default 0.6)
#   params[1] = fraction of the time with both reactors off (Default 0.0)
#   params[2] = fraction with which the rate in each bin gaussian fluctuates
def cns_time_simple(time_days,params):
    time_days = np.array(time_days)
    if(len(params)>0 and params[0]!=""):
        full_power_frac = params[0]
    else:
        full_power_frac = 0.6
    if(len(params)>1 and params[1]!=""):
        off_frac = params[1]
    else:
        off_frac = 0.0
    if(len(params)>2 and params[2]!=""):
        fluctuation_frac = params[2]
    else:
        fluctuation_frac = 0.05
    return cns_time_simple_0(time_days,full_power_frac,off_frac,fluctuation_frac)
# ----------------------------------------------------------------------

def cns_time_very_simple(time_days, frac=0.6):
    while time_days>365.0:
        time_days-=365.0
    if time_days<frac*365.0:
        return 1.0
    elif time_days<(frac+(1.0-frac)/2.0)*365.0:
        return 0.6
    else:
        return 0.4
