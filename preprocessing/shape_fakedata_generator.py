#!/usr/bin/evn python
#
# Preprocessing scripts to bin and store shapes and
# other info for a stan/morpho model
# ----------------------------------------------------------------------
# Author: J. P. Johnston <jpj13@mit.edu>
#
# The code to parse command line arguments and read data from .yaml
# files is based on morpho.py
# morpho is located at https://github.com/project8/morpho/tree/v1.1.5
# The authors of morpho are:
#   J. A. Formaggio <josephf@mit.edu>
#   T. E. Weiss <tweiss@mit.edu>
#   M. Guigue <mathieu.guigue@pnnl.gov>
#   J. N. Kofron <jared.kofron@gmail.com>
# ----------------------------------------------------------------------
#
# Creating .yaml files to work with shape_data_generator.py:
#
# renormalization: The options are "none", "integral", or "initial".
#   - "none" will simply rebin the given function or data without
#      renormalizing
#   - "integral" will renormalize so that the sum over all bins is 1.
#     This is useful if you know the total number of events over your
#     entire region of interest, such as with spectral shape
#   - "initial" will renormalize by dividing by the value of the
#     function at the lower bound. This is useful when your signal
#     magnitude is the initial signal, such as in the case of time
#     dependence.

import logging
logger = logging.getLogger('shapegenerator')
logger.setLevel(logging.DEBUG)
base_format = '%(asctime)s[%(levelname)-8s] %(name)s(%(lineno)d) -> %(message)s'
logging.basicConfig(format=base_format, datefmt='%m/%d/%Y %H:%M:%S')

from yaml import load as yload
from argparse import ArgumentParser
import pystan

import numpy as np
import ROOT as ROOT# import ROOT, TStyle, TCanvas, TH1F, TGraph, TLatex, TLegend, TFile, TTree, TGaxis, TRandom3, TNtuple, TTree
from array import array

import sys
from importlib import import_module

#import imp
#temp = imp.load_source('ric_shapefunctions','ric_functions/ric_shapefunctions.py')

# ----------------------------------------------------------------------
# Methods to parse command line arguments and read data from files
def parse_args():
    p = ArgumentParser(description='''
        Preprocessing of model parameters for a stan/morpho model
    ''')
    p.add_argument('-c','--config',
                   metavar='<configuration file>',
                   help='Full path to the configuration file',
                   required=True)
    p.add_argument('param',nargs='*',
                   default=False,
                   help='Manualy change of a parameter and its value')
    return p.parse_args()
def read_param(yaml_data, node, default):
        data = yaml_data
        xpath = node.split('.')
        try:
            for path in xpath:
                data = data[path]
        except Exception as exc:
            if default == 'required':
                err = """FATAL: Configuration parameter {0} required but not\
                provided in config file!
                """.format(node)
                logger.debug(err)
                raise exc
            else:
                data = default
        return data
# Input: root_file_path, tree_name, and branch_name
#        define where to read data
# Output: data is all data from the branches
def readTTree(root_file_path,tree_name,branch_name):
#    print('Reading {}'.format(root_file_path))
    myfile = ROOT.TFile(root_file_path,"READ")
    tree = myfile.Get(tree_name)
    n = int(tree.GetEntries())
    data = []
    for elt in tree:
        curr_data = getattr(elt,branch_name)
        data.append(curr_data)
    return data
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
               rescale_factor=1.0,tree_name="shape",
               x_name="xdata",y_name="ydata"):
    # Bin the shape
    dx = (ub-lb)/float(nbins)
    h = ROOT.TH1F("h","",nbins,lb,ub)
    hw = ROOT.TH1F("hw","",nbins,lb,ub)
    havg = ROOT.TH1F("havg","",nbins,lb,ub)
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
                      optional_outputs=[],label="signal_or_background"):
    # Settings about output messages and saved information
    print_debug = read_param(optional_outputs,'print_debug_statements',False)
    store_info = read_param(optional_outputs,'store_info_text',False)
    store_plots = read_param(optional_outputs,'store_info_plots',False)
    if(store_info or store_plots):
        info_dir = read_param(optional_outputs,'info_output_directory',"required")
    if(store_info):
        txt_file = open(info_dir + "/" + label + "_shape_info.txt","w")
        if(print_debug):
            print("Shape info will be stored in %s" % info_dir + "/" + label + "_shape_info.txt")
    if(print_debug and store_plots):
            print("Shape plots will be stored in %s" % info_dir)

    
    if(output_type != "root"):
        print("%s not yet implemented. Shapes currently are saved as \".root\".")
        
    # Create root file to save all shapes of this set of signals
    myfile = ROOT.TFile(output_file,"RECREATE")

    for i in range(0,len(vars_dict)):
        # iterate over all independent variables
        curr_var = vars_dict[i]
        curr_var_name = read_param(curr_var,'name','x_%i'%i)
        if(print_debug):
            print("Current independent variable: %s" % curr_var_name)
        nbins = read_param(curr_var,'bins',30)
        lb = read_param(curr_var,'lower_bound','required')
        ub = read_param(curr_var,'upper_bound','required')
        norm_type = read_param(curr_var,'renormalization','none')
        if norm_type == 'initial':
            norm_type = 'rescale' # Rescale by the initial value

        for idsig,sig in enumerate(signals_dict):
            # iterate over all signals
            sig_name = read_param(sig,'name','Signal_%i'%idsig)
            if(print_debug):
                print("\tCurrent %s: %s" % (label,sig_name))
            dim_params = read_param(sig,'dimension_params','required')
            dim_p = dim_params[i]

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
                
                # Call the function
                samples = read_param(dim_p,'samples_per_bin',10)
                x_vals = np.linspace(lb,ub,nbins*samples)
                fcn_params = read_param(dim_p,'params',[])
                if(len(fcn_params)>0):
                    y_vals = y_fcn(x_vals,fcn_params)
                else:
                    y_vals = y_fcn(x_vals)
                init_norm = y_fcn(lb)
            elif fcn_type == "data_files":
                data_file_type = read_param(dim_p,'data_file_type','csv')
                if(data_file_type!='csv' and data_file_type!='.csv'):
                    print("WARNING: data_file_type currently must be \"csv\"")
                    print("         data_file_type==\"%s\" invalid. Assuming data_file_type==\"csv\"." % data_file_type)
                x_path = read_param(dim_p,'x_data_location','required')
                y_path = read_param(dim_p,'y_data_location','required')
                x_vals = np.genfromtxt(x_path,delimiter=',')
                y_vals = np.genfromtxt(y_path,delimiter=',')
                init_norm = 1.0 # CHANGE THIS TO GET y(lb)
            elif fcn_type == "stan_fcn":
                x_vals = [1.0]
                y_vals = [1.0]
                init_norm = 1
                print("Reading shape from stan_fcn not yet implemented")
            else:
                print("ERROR: %s type==\"%s\" not valid" % (sig_name,fcn_type))
                print("       Type must be \"py_fcn\", \"data_files\", or \"stan_fcn\"")
                print("       Exiting shape generation for %s" % sig_name)
                return

            # Bin the x and y values into an nbin histogram
            tree_name = read_param(dim_p,'tree_name','required')
            curr_x_name = read_param(dim_p,'x_branch_name','xdata_%i'%i)
            curr_y_name = read_param(dim_p,'y_branch_name',sig_name)
            curr_tree = \
               rebin_data(x_vals,y_vals,nbins,lb,ub,\
                          norm_type,init_norm,\
                          tree_name,curr_x_name,curr_y_name)
            curr_tree.Write()

            # TODO: Implement storing debug text and plots
            # Print optional output text and plots about the current signal
            if(store_info):
                txt_file.write("To do: Implement storing info about the shape\n")
            if(store_plots):
                #curr_tree.Draw()
                #curr_tree.SaveAs(info_dir+"/"+out.pdf)
                if(print_debug):
                    print("\t\tTo do: Implement storing plots")
    # Close root and text files
    myfile.Close()
    if(store_info):
        txt_file.close()
    return
# ----------------------------------------------------------------------
# TODO: Allow gaussian distribution of each bin or globally for
#       each shape
def sum_weighted_shapes(input_root_file,shapes_dict,nBins,
                        curr_var_count,gauss_redist):
    # Get the sum of all signal shapes
    tot_shape = np.zeros(nBins)
    tot_weight = 0.0
    for sig in shapes_dict:
        curr_weight = read_param(sig,'fake_data_weight','required')
        tot_weight += curr_weight
        dim_p = read_param(sig,'dimension_params','required')
        curr_dim = dim_p[curr_var_count]
        tree_name = read_param(curr_dim,'tree_name','required')
        branch_name = read_param(curr_dim,'y_branch_name','required')
        curr_shape = np.array(readTTree(input_root_file,tree_name,branch_name))
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
#         The contents of append_R_file will be appended to the end of
#         output_file_name if output_file_type=="R". append_R_file
#         must be an R file. This allows multiple data arrays to be
# 	  stored in a single R file.
def write_data_array(arr,output_file_name,output_file_type,
                     tree_name,data_type="float",
                     append_R_file=""):
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
        print("ERROR: Writing to root file not yet implemented")
        if(ndim>1):
            print("ERROR: For >1 dimensional data, currently only storing to a \"R\" file is supported")
            print("       No data will be stored because %s is invalid" % output_file_type)
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
    if(store_info):
        txt_file = open(info_dir + "/fake_data_info.txt","w")
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
        # TODO: update total_shapes to read nBins from the additional_output_file
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

    ran = ROOT.TRandom3()
    ran.SetSeed(0)
    if(gauss_redist):
        # TODO: Take this input from the yaml file
        #       The 5% global flucutation accounts for us not knowing
        #       the overall power
        sig_mag = ran.Gaus(sig_mag,sig_mag*0.05)
        
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
            # TODO: Take this input from the .yaml file
            #       I use 15% because 15 eV is the resolution, bin
            #       width is 100 eV. I also need to think of a better
            #       way to address energy resolution, preferably one
            #       that maintains the total number of events
            curr_fake_data = ran.Gaus(curr_fake_data,0.15*curr_fake_data)
        if(poisson_redist):
            curr_fake_data = ran.Poisson(curr_fake_data)
        fake_data.itemset(index,curr_fake_data)
        
    # Write the array of fake data to file
    fake_data_output_file = read_param(fake_data_settings,'fake_data_output_file','required')
    fake_data_output_type = read_param(fake_data_settings,'fake_data_output_type','R')
    fake_data_output_tree = read_param(fake_data_settings,'fake_data_output_tree','fake_data')
    # If fake_data_output_type=="R", this will copy the contents of the
    # additional_stan_input file to the end of the fake_data file, so
    # only one R file is required by the stan model
    write_data_array(fake_data,fake_data_output_file,fake_data_output_type,
                     fake_data_output_tree,"int",
                     append_R_file=additional_stan_input_file)
    return
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

    vars_dict = read_param(prep_dict,'indep_vars','required')
    signal_dict = read_param(prep_dict,'signals','required')
    signal_output_file = read_param(prep_dict,'signal_shape_output_file','required')
    signal_output_type = read_param(prep_dict,'signal_shape_output_type',"root")
    back_dict = read_param(prep_dict,'backgrounds','required')
    back_output_file = read_param(prep_dict,'background_shape_output_file','required')
    back_output_type = read_param(prep_dict,'background_shape_output_type',"root")
    opt_out_settings = read_param(prep_dict,'optional_output_settings',[])

    # Generate and store the shapes of signals and backgrounds
    if read_param(prep_dict,'generate_shapes',True):
        print("Generating and storing signal and background shapes")
        create_shape_info(vars_dict,signal_dict,signal_output_file,signal_output_type,
                          opt_out_settings,"signal")
        create_shape_info(vars_dict,back_dict,back_output_file,back_output_type,
                          opt_out_settings,"background")
        # Store the number of bins so stan can access
        additional_file_name = read_param(prep_dict,'additional_stan_input_file','required')
        additional_file = open(additional_file_name,"w")
        for ivar,var in enumerate(vars_dict):
            nbins=read_param(var,'bins','required')
            lb = read_param(var,'lower_bound','required')
            ub = read_param(var,'upper_bound','required')
            additional_file.write('nBins_%i <- %i\n'%(ivar,nbins))
            additional_file.write('lb_%i <- %f\n'%(ivar,lb))
            additional_file.write('ub_%i <- %f\n'%(ivar,ub))
        additional_file.close()
        
    # Generate fake data
    if read_param(prep_dict,'generate_fake_data',False):
        print("Generating and storing fake data")
        additional_file_name = read_param(prep_dict,'additional_stan_input_file','required')
        generate_fake_data(vars_dict,signal_dict,signal_output_file,signal_output_type,
                           back_dict,back_output_file,back_output_type,
                           read_param(prep_dict,'fake_data_settings','required'),
                           opt_out_settings,additional_file_name)

    print("Preprocessing complete")
