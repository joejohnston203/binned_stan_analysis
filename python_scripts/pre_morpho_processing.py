#!/usr/bin/evn python
#
# Preprocessing scripts to bin and store shapes and
# other info for a stan/morpho model
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
#
# See scripts/example_ricochet_rate_analyzer.yaml or
# scripts/example_ricochet_rate_spectrum_analyzer.yaml as an example
# of how to format a .yaml file to work with this script.
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

from sensitivitytools.data_generation import *
from sensitivitytools.root_tools import *
from sensitivitytools.input_processing import *

# ----------------------------------------------------------------------
# Arguments:
#   xdata, ydata:  arrays of data points that will be binned
#   nbins, lb, ub: parameters for the returned tree of binned data
#   norm_type:     Options are 'none' (default), 'rescale', or 'integral'.
#                  'rescale' will divide by rescale_factor, and 'integral'
#                  will normalize so the sum over all bins is 1.
#   rescale_facotr: Factor used when norm_type='rescale'
#   tree_name,x_name,y_name: Properties for returned ttree
# Returns:
#   TTree with two arrays: One for the x values corresponding to the
#   bins, and one with the y values that gives the average value of the
#   function in each bin multiplied by the width of the bin. The y values
#   will be renormalized according to norm_type.
#
def rebin_data(xdata,ydata,nbins,lb,ub,norm_type="none",
               tree_name="shape",
               x_name="xdata",y_name="ydata",
               save_plot_name="",
               rescale_factor=1.0,
               variable_binning=None):

    # Bin the shape
    dx = (ub-lb)/float(nbins)
    if(variable_binning is None):
        h = ROOT.TH1F("h","",nbins,lb,ub)
        hw = ROOT.TH1F("hw","",nbins,lb,ub)
        havg = ROOT.TH1F("havg","",nbins,lb,ub)
    else:
        h = ROOT.TH1F("h","",len(variable_binning)-1,variable_binning)
        hw = ROOT.TH1F("hw","",len(variable_binning)-1,variable_binning)
        havg = ROOT.TH1F("havg","",len(variable_binning)-1,variable_binning)
    list_x = []
    list_y = []
    for j in range(0,len(xdata)):
        h.Fill(xdata[j],ydata[j]*dx)
        hw.Fill(xdata[j],1)
        # Bin numering starts at 1, not 0. List numbering starts at 0
    for j in range(1,nbins+1):
        list_x.append(h.GetBinCenter(j))
        list_y.append(h.GetBinContent(j)/max(1,hw.GetBinContent(j)))
        havg.Fill(h.GetBinCenter(j),list_y[j-1])
    # Properly normalize
    if(norm_type=="rescale"):
        # Normalize by dividing by the initial of the function
        list_y  = np.asarray(list_y)/float(rescale_factor)
    elif(norm_type=="integral"):
        # Normalize so the sum over all bins is 1
        sum_y = sum(list_y)
        list_y = np.asarray(list_y)/float(sum_y)
        # Otherwise do not renormalize, and assume that xdata and ydata
        # were already properly normalized
    # Create the TTree
    tmp_x = array('f',[ 0 ])
    tmp_y = array('f',[ 0 ])
    shape_tree = ROOT.TTree(tree_name,tree_name)
    shape_tree.Branch(x_name,tmp_x,x_name+'/F')
    shape_tree.Branch(y_name,tmp_y,y_name+'/F')
    # Write expected Data
    for j in range(0,len(list_x)):
        tmp_x[0] = list_x[j]
        tmp_y[0] = list_y[j]
        shape_tree.Fill()
    if(save_plot_name!=""):
        temp_gr = ROOT.TGraph(len(list_x),np.array(list_x),np.array(list_y))
        root_make_plot([temp_gr],save_plot_name,["APL"],tree_name,x_name,y_name,marker_styles=[3])
    return shape_tree
# ----------------------------------------------------------------------
# Method to accept a dictionary of signals or backgrounds, then
# generate and store shapes and info for the signals and backgrounds
# Inputs:
#  The first four inputs should be dictionaries or variables from a
#  .yaml file, as described at the top of this document
#  label: the type of shape, eg signal or background. This will be
#         used as a prefix for any optional output files that are saved,
#         and will be printed in debugging statements
def create_shape_info(vars_dict,signals_dict,output_file,output_type,
                      optional_outputs=[],label="signal_or_background",
                      additional_file_name=""):
    # open the additional file
    if(additional_file_name!=""):
        additional_file = open(additional_file_name,"a+")
    # Settings about output messages and saved information
    print_debug = read_param(optional_outputs,'print_debug_statements',False)
    store_info = read_param(optional_outputs,'store_info_text',False)
    store_plots = read_param(optional_outputs,'store_info_plots',False)
    if(store_info or store_plots):
        info_dir = read_param(optional_outputs,'info_output_directory',"required")
        create_path(info_dir,False)
    if(store_info):
        prefix = read_param(optional_outputs,'info_output_prefix','')
        if(prefix!=""):
            filename = info_dir + "/" + prefix+"_"+label + "_shape_info.txt"
        else:
            filename = info_dir + "/" + label + "_shape_info.txt"
        create_path(filename,True)
        txt_file = open(filename,"w")
        if(print_debug):
            print("%s info will be stored in %s" % (label,filename))
    if(print_debug and store_plots):
        print("%s plots will be stored in %s" % (label,info_dir))
    if(store_plots):
        tree_names = []
        x_branch_names = []
        y_branch_names = []

    
    if(output_type != "root"):
        print("%s not yet implemented. Shapes currently are saved as \".root\".")
    
    # Recreate root file used to save all shapes of this set of signals
    create_path(output_file,True)
    myfile = ROOT.TFile(output_file,"RECREATE")
    myfile.Close()
    if(store_info):
        txt_file.write("%s output file: %s\n"%(label,output_file))

    for var_index in range(0,len(vars_dict)):
        # iterate over all independent variables
        curr_var = vars_dict[var_index]
        curr_var_name = read_param(curr_var,'name','x_%i'%var_index)
        nbins = int(read_param(curr_var,'bins',30))
        lb = float(read_param(curr_var,'lower_bound','required'))
        ub = float(read_param(curr_var,'upper_bound','required'))
        binning_file = read_param(curr_var,'binning_file','none')
        if(binning_file!='none'):
            binning = np.loadtxt(binning_file)
        else:
            binning = None
        if(print_debug):
            print("Current independent variable: %s" % curr_var_name)
            print("\tlb = %f"%lb)
            print("\tub = %f"%ub)
            print("\tnBins = %i"%nbins)
        if(store_info):
            txt_file.write("Current independent variable: %s\n" % curr_var_name)
            txt_file.write("\tlb = %f\n"%lb)
            txt_file.write("\tub = %f\n"%ub)
        norm_type = read_param(curr_var,'renormalization','none')
        if norm_type == 'initial':
            norm_type = 'rescale' # Rescale by the initial value

        for idsig,sig in enumerate(signals_dict):
            # iterate over all signals
            sig_name = read_param(sig,'name','Signal_%i'%idsig)
            if(print_debug):
                print("\tCurrent %s: %s" % (label,sig_name))
            dim_params = read_param(sig,'dimension_params','required')
            dim_p = dim_params[var_index]

            fcn_type = read_param(dim_p,'type','required')

            # Get x and y values for the current function
            if fcn_type == "py_fcn":
                # Get the function
                path=read_param(dim_p,'location','required')
                mod_name=read_param(dim_p,'module','required')
                fcn_name=read_param(dim_p,'fcn_name','required')
                sys.path.append(path)
                temp_module = import_module(mod_name)
                y_fcn=getattr(temp_module,fcn_name)
                
                # Define x_values so it is the center value of each bin
                samples = int(read_param(dim_p,'samples_per_bin',10))
                bin_width = (ub-lb)/float(nbins)
                x_vals = np.linspace(lb,ub,nbins*samples,False)+bin_width/2.0
                fcn_params = read_param(dim_p,'params',[])
                if(len(fcn_params)>0):
                    y_vals = y_fcn(x_vals,fcn_params)
                else:
                    y_vals = y_fcn(x_vals)
            elif fcn_type == "data_files":
                data_file_type = read_param(dim_p,'data_file_type','csv')
                if(data_file_type=='csv'):
                    x_path = read_param(dim_p,'x_data_location','required')
                    x_vals = np.genfromtxt(x_path,delimiter=',')
                    y_path = read_param(dim_p,'y_data_location','required')
                    y_vals = np.genfromtxt(y_path,delimiter=',')
                elif(data_file_type=='columns'):
                    path = read_param(dim_p,'data_location','required')
                    delim = read_param(dim_p,'delimiter','\t')
                    data = np.loadtxt(path,delimiter=delim)
                    x_vals = data[:,0]
                    y_vals = data[:,1]
                else:
                    print("WARNING: data_file_type currently must be \"csv\" or \"columns\"")
                    print("         data_file_type==\"%s\" invalid. Exiting." % data_file_type)
                    return
            elif fcn_type == "stan_fcn":
                x_vals = [1.0]
                y_vals = [1.0]
                print("Reading shape from stan_fcn not yet implemented")
            elif fcn_type == "root_histo":
                tempfile = ROOT.TFile(read_param(dim_p, 'location', 'required'), "READ")
                tree = tempfile.Get(read_param(dim_p, 'histo_tree', 'required'))
                if(binning is None):
                    print("Fixed binning")
                    rh = ROOT.TH1F("rh","",nbins,lb,ub)
                else:
                    print("Variable binning")
                    rh = ROOT.TH1F("rh","",len(binning)-1,binning)
                br_name = read_param(dim_p, 'histo_branch', 'required')
                tree.Draw(br_name+">>rh","","goff")
                x_vals = np.empty(len(binning)-1)
                y_vals = np.empty(len(binning)-1)
                for i in range(0,len(binning)-1):
                    # np arr index starts at 0, root index starts at 1
                    x_vals[i] = rh.GetBinCenter(i+1);
                    y_vals[i] = rh.GetBinContent(i+1);
                tempfile.Close()
            else:
                print("ERROR: %s type==\"%s\" not valid" % (sig_name,fcn_type))
                print("       Type must be \"py_fcn\", \"data_files\", or \"stan_fcn\"")
                print("       Exiting shape generation for %s" % sig_name)
                return

            # Bin the x and y values into an nbin histogram
            tree_name = read_param(dim_p,'tree_name','required')
            if(store_plots):
                save_plot_file_name = info_dir+"/"+tree_name+".pdf"
            else:
                save_plot_file_name = ""
            curr_x_name = read_param(dim_p,'x_branch_name','xdata_%i'%i)
            curr_y_name = read_param(dim_p,'y_branch_name',sig_name)

            curr_tree = \
               rebin_data(x_vals,y_vals,len(x_vals),lb,ub,\
                          norm_type,\
                          tree_name,curr_x_name,curr_y_name,
                          save_plot_file_name,
                          variable_binning=binning)
            myfile = ROOT.TFile(output_file,"UPDATE")
            curr_tree.Write()
            myfile.Close()

            # Access and store the fractions for gaussian redistribution
            if(additional_file_name!=""):
                curr_bins_frac = read_param(dim_p,'bin_gauss_var_frac',0.0)
                curr_bins_var = read_param(dim_p,'stored_bin_frac_name',"")
                if(curr_bins_var!=""):
                    additional_file.write(curr_bins_var+" <- %f\n"%curr_bins_frac)
                curr_global_frac = read_param(dim_p,'global_gauss_var_frac',0.0)
                curr_global_var = read_param(dim_p,'stored_global_frac_name',"")
                if(curr_global_var!=""):
                    additional_file.write(curr_global_var+" <- %f\n"%curr_global_frac)

            
            if(store_info):
                txt_file.write("--Current %s: %s\n" % (label,sig_name))
                txt_file.write("\tFunction type = %s\n" % fcn_type)
                txt_file.write("\ttree name = %s\n" % tree_name)
                txt_file.write("\tx branch name = %s\n" % curr_x_name)
                txt_file.write("\ty branch name = %s\n" % curr_y_name)
                if(curr_bins_var!=""):
                    txt_file.write("\tfake data bins gauss frac = %.3e, variable = %s\n" % (curr_bins_frac,curr_bins_var))
                if(curr_global_var!=""):
                    txt_file.write("\tfake data global gauss frac = %.3e, variable = %s\n" % (curr_global_frac,curr_global_var))

    if(store_info):
        txt_file.close()
    if(additional_file_name!=""):
        additional_file.close()
    return
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
# Note that I read all text from both files into memory, so this
# method may crash for very large files
def concatenate_files(file1_name, file2_name,output_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    text = file1.read() + "\n" + file2.read()
    file1.close()
    file2.close()
    outfile = open(output_name,'w')
    outfile.write(text)
    outfile.close()
    
# ----------------------------------------------------------------------
# Inputs: numpy array of data, path to output file, and output type.
#         Currently, the output type can only be "R" or "root", and
#         for an array of dimension>=2 only "R" is supported
#
#         branch_name is only used for root output files.
#
#         The contents of append_R_file will be appended to the end of
#         output_file_name if output_file_type=="R". append_R_file
#         must be an R file. This allows multiple data arrays to be
# 	  stored in a single R file.
def write_data_array(arr,output_file_name,output_file_type,
                     tree_name,branch_name="data",
                     data_type="float",
                     append_R_file=""):
    create_path(output_file_name,True)
    ndims = len(arr.shape)
    if(output_file_type=='R'):
        # Store data to an R file
        if(data_type=="float"):
            tempdict = {tree_name:arr.astype(float)}
        elif(data_type=="int"):
            tempdict = {tree_name:arr.astype(int)}
        else:
            print("ERROR: Invalid data_type %s- select \"float\" or \"int\"" % data_type)
        pystan.misc.stan_rdump(tempdict, output_file_name)
        # Append the contents of append_R_file to the output_file
        if(append_R_file != ""):
            concatenate_files(output_file_name,append_R_file,output_file_name)
    elif(output_file_type=='root'):
        if(ndims>1):
            print("ERROR: For >1 dimensional data, currently only storing to a \"R\" file is supported")
            print("       No data will be stored because %s is invalid" % output_file_type)
# vars: output_file_name, tree_name, branch_name
        myfile = ROOT.TFile(output_file_name,"RECREATE")
        curr_tree = ROOT.TTree(tree_name,tree_name)
        if(data_type=="int"):
            tmp_arr = array('i',[ 0 ])
            curr_tree.Branch(branch_name,tmp_arr,branch_name+'/I')
        else:
            tmp_arr = array('f',[ 0 ])
            curr_tree.Branch(branch_name,tmp_arr,branch_name+'/F')
        for i in range(0,len(arr)):
            tmp_arr[0] = arr[i]
            curr_tree.Fill()
        curr_tree.Write()
        myfile.Close()
    else:    
        print("ERROR: %s is not a valid fake data output type. Select \"R\" or \"root\". No data will be stored."
              % output_file_type)
    return
# ----------------------------------------------------------------------
# Method to generate fake data
# Input: Dictionaries from a yaml file as described at the top of this
#        file. This method uses that info to find the shapes stored by
#        create_shape_info()
# stored 
# Output: Creates an n-dimensional array of fake data. This is then
#         stored, currently as a "root" file for one dimensional data,
#         and in an "R" file for two dimensional data. Higher
#         dimensional data currently cannot be stored.
def generate_fake_data(vars_dict,signals_dict,signals_input_file,signals_input_type,
                       backs_dict,backs_input_file,backs_input_type,
                       fake_data_settings,optional_outputs=[],
                       additional_stan_input_file=""):
    # Settings about output messages and saved information
    print_debug = read_param(optional_outputs,'print_debug_statements',False)
    store_info = read_param(optional_outputs,'store_info_text',False)
    store_plots = read_param(optional_outputs,'store_info_plots',False)
    if(store_info or store_plots):
        info_dir = read_param(optional_outputs,'info_output_directory',"required")
        create_path(info_dir,False)

    if(store_info):
        prefix = read_param(optional_outputs,'info_output_prefix','')
        if(prefix!=""):
            filename = info_dir + "/" + prefix+"_fake_data_info.txt"
        else:
            filename = info_dir + "/fake_data_info.txt"
        create_path(filename,True)
        txt_file = open(filename,"w")
        if(print_debug):
            print("Fake data info will be stored in %s" % info_dir + "/fake_data_info.txt")
    if(print_debug and store_plots):
            print("Fake data plots will be stored in %s" % info_dir)
    
    if(signals_input_type != "root" or backs_input_type != "root"):
        print("%s not yet implemented for shapes. Shapes currently are saved as \".root\".")

    gauss_redist = read_param(fake_data_settings,'fake_gaussian_redistribution',True)
    poisson_redist = read_param(fake_data_settings,'fake_poisson_redistribution',True)

    
    # Create arrays of total shape for each independent variable
    total_sig_shapes = list()
    total_back_shapes = list()
    nBins_arr = np.empty(len(vars_dict),dtype=int)
    for i in range(0,len(vars_dict)):
        curr_var=vars_dict[i]
        curr_var_name = read_param(curr_var,'name','x_%i'%i)
        if(print_debug):
            print("Current independent variable: %s" % curr_var_name)
        nBins = read_param(curr_var,'bins','required')
        nBins_arr[i] = np.int(nBins)
        total_sig_shapes.append(sum_weighted_shapes(signals_input_file,signals_dict,nBins,i,gauss_redist))
        total_back_shapes.append(sum_weighted_shapes(backs_input_file,backs_dict,nBins,i,gauss_redist))
        
    # Generate an n-dimensional array of fake data
    sig_mag = read_param(fake_data_settings,'fake_signal_magnitude','required')
    back_mag = read_param(fake_data_settings,'fake_background_magnitude','required')
    if(print_debug):
        print("Fake data signal magnitude = %f"%sig_mag)
        print("Fake data background magnitude = %f"%back_mag)

    fake_data = np.empty(nBins_arr)
    ran = ROOT.TRandom3()
    ran.SetSeed(0)
    for index,val in np.ndenumerate(fake_data):
        curr_sig_counts = sig_mag
        curr_back_counts = back_mag
        for j in range(0,len(vars_dict)):
            curr_sig_counts*=total_sig_shapes[j][index[j]]
            curr_back_counts*=total_back_shapes[j][index[j]]
        curr_fake_data = curr_sig_counts+curr_back_counts
        if(gauss_redist):
            gauss_frac = read_param(fake_data_settings,'bin_gauss_var_total',0.0)
            if(gauss_frac!=0.0):
                curr_fake_data = ran.Gaus(curr_fake_data,gauss_frac*curr_fake_data)
        if(poisson_redist):
            curr_fake_data = ran.Poisson(curr_fake_data)
        fake_data.itemset(index,curr_fake_data)
        
    # Write the array of fake data to file
    fake_data_output_file = read_param(fake_data_settings,'fake_data_output_file','required')
    fake_data_output_type = read_param(fake_data_settings,'fake_data_output_type','R')
    fake_data_output_tree = read_param(fake_data_settings,'fake_data_output_tree','fake_data')
    fake_data_output_branch = read_param(fake_data_settings,'fake_data_output_branch','fake_data_counts')
    # If fake_data_output_type=="R", this will copy the contents of the
    # additional_stan_input file to the end of the fake_data file, so
    # only one R file is required by the stan model
    write_data_array(fake_data,fake_data_output_file,fake_data_output_type,
                     fake_data_output_tree,data_type="int",
                     append_R_file=additional_stan_input_file)
    if(store_info):
        txt_file.write("Total Signal Magnitude = %.3e\n"%sig_mag)
        txt_file.write("Total Background Magnitude = %.3e\n"%back_mag)
        txt_file.write("Gaussian Redistribution = %s\n"%gauss_redist)
        # TODO: write gaussian redist fractions
        txt_file.write("Poisson Redistribution = %s\n"%poisson_redist)
        txt_file.close()
    if(store_plots):
        # Get the x values for each dimension
        x_file = signals_input_file
        sig = signals_dict[0]
        dim_params = read_param(sig,'dimension_params','required')
        for i,dim_p in enumerate(dim_params):
            x_tree = read_param(dim_p,'tree_name','required')
            x_branch = read_param(dim_p,'x_branch_name','required')
            x_data = readTTree(x_file,x_tree,x_branch)
            if(len(fake_data.shape)==1):
                y_data = fake_data
            elif(len(fake_data.shape)==2):
                if(i==0):
                    #y_data = fake_data[:,0]
                    y_data = np.sum(fake_data,1)
                else:
                    #y_data = fake_data[0,:]
                    y_data = np.sum(fake_data,0)
            elif(len(fake_data.shape)>2):
                print("Fake data plotting not implemented for more than 2 dimensional data")
            plot_file_name = info_dir + "/fake_data_%i.pdf"%i
            
            temp_gr = ROOT.TGraph(len(x_data),np.array(x_data),np.array(y_data))
            root_make_plot([temp_gr],plot_file_name,["APL"],"Fake Data",x_branch,"Counts",marker_styles=[3])
            
    return
# ----------------------------------------------------------------------
# Main method
# Input: preprocessing dictionary
# Can generate and store shapes for each signal and background
# Can generate fake data
def process(prep_dict):
    vars_dict = read_param(prep_dict,'indep_vars','required')
    signal_dict = read_param(prep_dict,'signals','required')
    signal_output_file = read_param(prep_dict,'signal_shape_output_file','required')
    signal_output_type = read_param(prep_dict,'signal_shape_output_type',"root")
    back_dict = read_param(prep_dict,'backgrounds','required')
    back_output_file = read_param(prep_dict,'background_shape_output_file','required')
    back_output_type = read_param(prep_dict,'background_shape_output_type',"root")
    opt_out_settings = read_param(prep_dict,'optional_output_settings',[])

    set_root_env() # Set ROOT to make nice plots, in case any plots are saved

    # Generate and store the shapes of signals and backgrounds
    if read_param(prep_dict,'generate_shapes',True):
        print("Generating and storing signal and background shapes")
        
        # Store the number of bins so stan can access, recreating
        # the additional file in the process
        additional_file_name = read_param(prep_dict,'additional_stan_input_file','required')
        create_path(additional_file_name,True)
        additional_file = open(additional_file_name,"w")
        for ivar,var in enumerate(vars_dict):
            nbins=read_param(var,'bins','required')
            lb = read_param(var,'lower_bound','required')
            ub = read_param(var,'upper_bound','required')
            additional_file.write('nBins_%i <- %i\n'%(ivar,nbins))
            additional_file.write('lb_%i <- %f\n'%(ivar,lb))
            additional_file.write('ub_%i <- %f\n'%(ivar,ub))
        fake_data_settings=read_param(prep_dict,'fake_data_settings',[])
        total_gauss_frac_var_name = read_param(fake_data_settings,'bin_gauss_var_name_total',"")
        total_gauss_frac_bins = read_param(fake_data_settings,'bin_gauss_var_total',0.0)
        if(total_gauss_frac_var_name!=""):
            additional_file.write(total_gauss_frac_var_name+" <- %f\n"%total_gauss_frac_bins)
        additional_file.close()

        create_shape_info(vars_dict,signal_dict,signal_output_file,signal_output_type,
                          opt_out_settings,"signal",additional_file_name)
        create_shape_info(vars_dict,back_dict,back_output_file,back_output_type,
                          opt_out_settings,"background",additional_file_name)
        
    # Generate fake data
    if read_param(prep_dict,'generate_fake_data',False):
        print("Generating and storing fake data")
        additional_file_name = read_param(prep_dict,'additional_stan_input_file','required')
        generate_fake_data(vars_dict,signal_dict,signal_output_file,signal_output_type,
                           back_dict,back_output_file,back_output_type,
                           read_param(prep_dict,'fake_data_settings','required'),
                           opt_out_settings,additional_file_name)

    print("Preprocessing complete")
# ----------------------------------------------------------------------
# Main method
if __name__== '__main__':
    # Get settings from .yaml configuration file
    args = parse_args()
    with open(args.config, 'r') as cfile:
        try:
            cdata = yload(cfile)
            if args.param:
                cdata = update_from_arguments(cdata,args.param)
        except Exception as err:
            logger.debug(err)
    prep_dict = read_param(cdata,'preprocessing','required')
    process(prep_dict)
# ----------------------------------------------------------------------
