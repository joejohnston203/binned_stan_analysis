# !/usr/bin/env python
#
# Code to take an input .yaml configuration file, run a frequentist
# sensitivity analysis, and create plots
# ----------------------------------------------------------------------

'''import logging
logger = logging.getLogger('shapegenerator')
logger.setLevel(logging.DEBUG)
base_format = '%(asctime)s[%(levelname)-8s] %(name)s(%(lineno)d) -> %(message)s'
logging.basicConfig(format=base_format, datefmt='%m/%d/%Y %H:%M:%S')

import pystan

import numpy as np
import ROOT as ROOT# import ROOT, TStyle, TCanvas, TH1F, TGraph, TLatex, TLegend, TFile, TTree, TGaxis, TRandom3, TNtuple, TTree
from array import array

import sys
from importlib import import_module
import os

from pre_morpho_processing import *
from sensitivitytools.root_tools import *
from sensitivitytools.input_processing import *'''

import numpy as np

#import sensitivitytools.frequentist_analysis as fa
from sensitivitytools.frequentist_analysis import *

if __name__=='__main__':
    dim_time = Dimension(name="time", lb=0, ub=365, nbins=5,
                         multiply_rate=True, convolution_fcn=None)
    dim_energy = Dimension(name="energy",lb=0,ub=1000,nbins=5,
                           multiply_rate=False, convolution_fcn=None)

    def f1(x):
        year = 365
        while(x>year):
            x -= 365
        if(x<0.6*year):
            return 1.0
        else:
            return 0.5
    f1 = np.vectorize(f1)
    shf1 = Shape(name="cns time", fcn=f1, fcn_args=[],
                 dimension_name="time", samples_per_bin=100)

    def f2(x):
        return np.exp(-x/1000.0)
    f2 = np.vectorize(f2)
    shf2 = Shape(name="cns energy", fcn=f2, fcn_args=[],
                 dimension_name="energy", samples_per_bin=100)

    parf = Parameter(name="cns_rate", shapes=[shf1,shf2], renormalize=True,
                     nuisance=False,
                     gauss_constraint=False, gauss_mean=5.0, gauss_sigma=0.25,
                     magnitude=5.0, magnitude_H0=0.0, lb=0.0, ub=50.0)

    def g1(x):
        return 1.0
    g1 = np.vectorize(g1)
    shg1 = Shape(name="back time", fcn=g1, fcn_args=[],
                 dimension_name="time", samples_per_bin=100)


    def g2(x):
        return 1.0
    g2 = np.vectorize(g2)
    shg2 = Shape(name="back energy", fcn=g2, fcn_args=[],
                 dimension_name="energy", samples_per_bin=100)

    parg = Parameter(name="back_rate", shapes=[shg1,shg2], renormalize=True,
                     nuisance=True,
                     gauss_constraint=False, gauss_mean=5.0, gauss_sigma=0.25,
                     magnitude=14.3, magnitude_H0=None, lb=0.0, ub=40.0)
    params = [parf, parg]

    def tot_rate(params, expt_settings=None, dims=None):
        return 0.5*np.sum(params)

    my_dims = [dim_time,dim_energy]

    expt1 = Experiment(name="Zn cns expt", dimensions=my_dims,
                       parameters=params, binned_expt=True,
                       expt_settings=None,
                       tot_rate_fcn=tot_rate, data=None)


    expt2 = Experiment(name="Zn cns another expt", dimensions=my_dims,
                       parameters=params, binned_expt=True,
                       expt_settings=None,
                       tot_rate_fcn=tot_rate, data=None)

    expts = [expt1,expt2]

    fa = FrequentistAnalysis("my analysis", binned_analysis=True,
                             experiments=expts, parameters=params)


    '''fa.frequentist_analysis("fa.root",200)
    print(fa)

    fa.save_all_shapes_pdf("outputs/shapes")
    fa.save_all_shapes_root("outputs/shapes")
    fa.save_all_data_root("outputs/data.root")
    fa.save_all_data_r("outputs/r_data")
    fa.save_all_data_pdf("outputs/data_pdfs")'''

    fa.generate_all_shapes()
    fa.generate_fake_data()
    print(fa.maximize_log_likelihood(use_H0=True, profile_nuisance=True))
    #print(fa.maximize_log_likelihood(use_H0=False))
