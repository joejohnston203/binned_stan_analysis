#!/usr/bin/evn python
#
# Methods to bin and store shapes and
# other info, as well as generate fake data for a sensitivity analysis
# ----------------------------------------------------------------------
# Author: J. P. Johnston <jpj13@mit.edu>
#
# The code to parse command line arguments and read data from .yaml
# files is from morpho.py
# morpho is located at https://github.com/project8/morpho/tree/v1.1.5
# The authors of morpho are:
#   J. A. Formaggio <josephf@mit.edu>
#   T. E. Weiss <tweiss@mit.edu>
#   M. Guigue <mathieu.guigue@pnnl.gov>
#   J. N. Kofron <jared.kofron@gmail.com>
# ----------------------------------------------------------------------

import logging
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

from .root_tools import *
from .input_processing import *

# ----------------------------------------------------------------------
# Takes inputs from a yaml file that indicates the locations of each
# signal or background shape inside the input root file, and computes
# a weighted average of the shapes.
#
# If weights is not empty, it is used to calculate the weighted
# average instead of the values in shapes_dict
def sum_weighted_shapes(input_root_file,shapes_dict,nBins,
                        curr_var_count,gauss_redist,
                        weights=[]):
    # Get the sum of all signal shapes
    tot_shape = np.zeros(nBins)
    tot_weight = 0.0
    for index, sig in enumerate(shapes_dict):
        if(index<len(weights)):
            curr_weight = weights[index]
        else:
            curr_weight = read_param(sig,'fake_data_weight','required')
        tot_weight += curr_weight
        dim_p = read_param(sig,'dimension_params','required')
        curr_dim = dim_p[curr_var_count]
        tree_name = read_param(curr_dim,'tree_name','required')
        branch_name = read_param(curr_dim,'y_branch_name','required')
        curr_shape = np.array(readTTree(input_root_file,tree_name,branch_name))
        if(gauss_redist):
            ran = ROOT.TRandom3()
            ran.SetSeed(0)
            curr_bins_frac = read_param(curr_dim,'bin_gauss_var_frac',0.0)
            if(curr_bins_frac!=0.0):
                for i in range(0,len(curr_shape)):
                    curr_shape[i] = ran.Gaus(curr_shape[i],curr_bins_frac*curr_shape[i])
            curr_global_frac = read_param(curr_dim,'global_gauss_var_frac',0.0)
            if(curr_global_frac!=0.0):
                curr_weight = ran.Gaus(curr_weight,curr_weight*curr_global_frac)
        tot_shape = tot_shape + curr_weight*curr_shape
    tot_shape = tot_shape/tot_weight
    return tot_shape
# ----------------------------------------------------------------------
