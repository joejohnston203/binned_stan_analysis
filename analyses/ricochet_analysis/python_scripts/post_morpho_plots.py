#!/usr/bin/evn python
#
# Plotting script to create histograms and other plots from the
# outputs of a stan/morpho run
# ----------------------------------------------------------------------
# Author: J. P. Johnston <jpj13@mit.edu>
#
# ----------------------------------------------------------------------
#
# See scripts/example_ricochet_rate_analyzer.yaml or
# scripts/example_ricochet_rate_spectrum_analyzer.yaml as an example of
# how to format a .yaml file to work with this script. Note that
# pre_morpho_processing.py must be run before some of these plots
# can be made.
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

import sensitivitytools.root_tools
from sensitivitytools.root_tools import *
from sensitivitytools.input_processing import *

# BEFORE COMITTING: Remove pre_morpho_processsing
from sensitivitytools.data_generation import sum_weighted_shapes

# ----------------------------------------------------------------------
def make_hist_plot(hist_settings,out_dir,out_prefix="",title_postfix=""):
    out_path = out_dir + out_prefix +\
                 read_param(hist_settings,'output_name','hist')
    title = read_param(hist_settings,'plot_title','')+title_postfix
    x_label = read_param(hist_settings,'x_label','')
    y_label = read_param(hist_settings,'y_label','')
    branches = read_param(hist_settings,'branches','required')
    file_names = []
    tree_names = []
    branch_names = []
    labels = []
    colors = [ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kMagenta,
              ROOT.kCyan,ROOT.kOrange,ROOT.kGreen+4,ROOT.kBlack]
    for i,b in enumerate(branches):
        file_names.append(read_param(b,'root_file_name','required'))
        tree_names.append(read_param(b,'tree_name','required'))
        branch_names.append(read_param(b,'branch_name','required'))
        labels.append(read_param(b,'label',""))
        curr_color = read_param(b,'color',"")
        if i>=len(colors):
            colors.append(ROOT.kBlack)
        if curr_color!="":
            colors[i] = curr_color
    add_legend = False
    for l in labels:
        if l!="":
            add_legend = True
    if(not add_legend):
        labels = []
    lxs = read_param(hist_settings,'leg_xstart',0.7)
    lys = read_param(hist_settings,'leg_ystart',0.75)
    lxe = read_param(hist_settings,'leg_xend',0.99)
    lye = read_param(hist_settings,'leg_yend',0.9)

    quantiles = read_param(hist_settings,'print_quantile_fracs',[])

    print_mean_stddev = read_param(hist_settings,'print_mean_stddev',True)
    root_plot_histogram(file_names,tree_names,branch_names,
                        out_path,title,x_label,y_label,
                        labels,lxs,lys,lxe,lye,
                        colors,quantiles,print_mean_stddev)
    return
# ----------------------------------------------------------------------
def make_corr_plot(corr_settings,out_dir,out_prefix="",title_postfix=""):
    file_name_x = read_param(corr_settings,'root_file_name_x','required')
    tree_name_x = read_param(corr_settings,'tree_name_x','required')
    branch_name_x = read_param(corr_settings,'branch_name_x','required')
    label_x = read_param(corr_settings,'x_label',"")
    file_name_y = read_param(corr_settings,'root_file_name_y','required')
    tree_name_y = read_param(corr_settings,'tree_name_y','required')
    branch_name_y = read_param(corr_settings,'branch_name_y','required')
    label_y = read_param(corr_settings,'y_label',"")
    title = read_param(corr_settings,'plot_title',"") + title_postfix

    out_path = out_dir + "/" + out_prefix +\
                read_param(corr_settings,'output_name','corr.pdf')
    root_make_corr_plot(file_name_x,tree_name_x,branch_name_x,
                        file_name_y,tree_name_y,branch_name_y,
                        out_path,label_x,label_y,title)
    return
# ----------------------------------------------------------------------
def extract_means(root_file_name,tree_names,branch_names):
    means=[]
    if(len(tree_names)!=len(branch_names)):
        return means

    for i in range(0,len(tree_names)):
        curr_file = ROOT.TFile(root_file_name)
        curr_tree = curr_file.Get(tree_names[i])
        curr_tree.Draw(branch_names[i]+">>temp_hist","","goff")
        temp_h = ROOT.TH1F(ROOT.gDirectory.Get("temp_hist"))
        means.append(temp_h.GetMean())
        curr_file.Close()
    return means
# ----------------------------------------------------------------------
def make_data_plot(data_settings,out_dir,out_prefix="",
                               title_postfix="",prep_dict=[]):
    out_path = out_dir + "/" + out_prefix +\
               read_param(data_settings,'output_name','data.pdf')
    title = read_param(data_settings,'plot_title',"")+title_postfix
    label_x = read_param(data_settings,'x_label',"")
    label_y = read_param(data_settings,'y_label',"")
    lxs = read_param(data_settings,'leg_xstart',0.7)
    lys = read_param(data_settings,'leg_ystart',0.75)
    lxe = read_param(data_settings,'leg_xend',0.99)
    lye = read_param(data_settings,'leg_yend',0.9)

    # Get x data and store it in a TChain object
    x_file_name = read_param(data_settings,'x_points_file_name','required')
    x_type = read_param(data_settings,'x_points_type','required')
    x_tree = read_param(data_settings,'x_points_tree_name','required')
    x_branch = read_param(data_settings,'x_points_branch_name','required')

    if x_type=="root":
        x_arr = readTTree(x_file_name,x_tree,x_branch)
    elif x_type=="R":
        x_dict = pystan.misc.read_rdump(x_file_name)
        x_arr = read_param(x_dict,x_tree,'required')
        x_elts = read_param(x_elts,'x_elts',[])
        if(len(x_elts)>0):
            x_arr = return_elts_of_arr(x_arr,x_elts)
    else:
        print("%s is not a valid type, please select \"root\" or \"R\""%x_type)
    x_arr = np.array(x_arr)
    n = len(x_arr)
    
    ymin = read_param(data_settings,'y_min',"")
    ymax = read_param(data_settings,'y_max',"")

    colors=[]
    marker_styles=[]
    draw_opts=[]
    leg_opts=[]
    labels=[]
    objs_to_print=[]
    text_to_print=""
    curr_var_num = read_param(data_settings,'curr_var',0)
    for (i,y_dict) in enumerate(read_param(data_settings,'y_points','required')):
        colors.append(read_param(y_dict,'color',ROOT.kBlack))
        marker_styles.append(read_param(y_dict,'marker',ROOT.kFullCircle))
        draw_opts.append(read_param(y_dict,'draw_opt',""))
        if i==0 and not 'A' in draw_opts[i]:
            draw_opts[i] = 'A' + draw_opts[i]
        leg_opts.append(read_param(y_dict,'leg_opt',""))
        labels.append(read_param(y_dict,'label',""))

        computed_points = read_param(y_dict,'computed_points',"")
        err_bar_type = read_param(y_dict,'error_bar_type',"none")
        err_bar_val = read_param(y_dict,'error_bar_val',0)
        y_arr = np.zeros(n)
        y_err = np.zeros(n)
        if(computed_points=="true_data" or computed_points=="extracted_data"):
            vars_dict = read_param(prep_dict,'indep_vars','required')
            curr_var = vars_dict[curr_var_num]
            signal_output_file = read_param(prep_dict,'signal_shape_output_file','required')
            signals_dict = read_param(prep_dict,'signals','required')
            back_output_file = read_param(prep_dict,'background_shape_output_file','required')
            back_dict = read_param(prep_dict,'backgrounds','required')
            if(computed_points=="true_data"):
                tot_sig_shape = sum_weighted_shapes(signal_output_file,signals_dict,n,curr_var_num,False)
                tot_back_shape = sum_weighted_shapes(back_output_file,back_dict,n,curr_var_num,False)
                fake_data_dict = read_param(prep_dict,'fake_data_settings','required')
                sig_mag = read_param(fake_data_dict,'fake_signal_magnitude','required')
                back_mag = read_param(fake_data_dict,'fake_background_magnitude','required')
                renorm=1.0
                if(len(vars_dict)==2):
                    # The displayed fake data will be summed in the
                    # other dimension, so we need to rescale
                    other_var_num = (curr_var_num+1)%2
                    other_var = vars_dict[other_var_num]
                    other_n = read_param(other_var,'bins',30)
                    other_sig_shape = sum_weighted_shapes(signal_output_file,signals_dict,other_n,other_var_num,False)
                    other_back_shape = sum_weighted_shapes(back_output_file,back_dict,other_n,other_var_num,False)
                    renorm = (sig_mag*np.sum(other_sig_shape)+back_mag*np.sum(other_back_shape))/(sig_mag+back_mag)
                y_arr = sig_mag*tot_sig_shape+back_mag*tot_back_shape
                y_arr = renorm*(sig_mag*tot_sig_shape+back_mag*tot_back_shape)
            else:
                # Extract signal weights
                sig_file_name = read_param(y_dict,'signal_distribution_file','required')
                sig_weight_trees = read_param(y_dict,'signal_weight_trees',[])
                sig_weight_branches = read_param(y_dict,'signal_weight_branches',[])
                sig_weights = extract_means(sig_file_name,sig_weight_trees,sig_weight_branches)
                tot_sig_shape = sum_weighted_shapes(signal_output_file,signals_dict,n,curr_var_num,False,sig_weights)
                # Extract signal magnitude
                sig_tree_name = read_param(y_dict,'signal_distribution_tree','required')
                sig_branch_name = read_param(y_dict,'signal_distribution_branch','required')
                sig_file = ROOT.TFile(sig_file_name)
                sig_tree = sig_file.Get(sig_tree_name)
                sig_tree.Draw(sig_branch_name+">>temp_hist","","goff")
                sig_hist = ROOT.TH1F(ROOT.gDirectory.Get("temp_hist"))
                sig_mag = sig_hist.GetMean()

                # Extract background weights
                back_file_name = read_param(y_dict,'background_distribution_file','required')
                back_weight_trees = read_param(y_dict,'background_weight_trees',[])
                back_weight_branches = read_param(y_dict,'background_weight_branches',[])
                back_weights = extract_means(back_file_name,back_weight_trees,back_weight_branches)
                tot_back_shape = sum_weighted_shapes(back_output_file,back_dict,n,curr_var_num,False,back_weights)
                # Extract background magnitude
                back_tree_name = read_param(y_dict,'background_distribution_tree','required')
                back_branch_name = read_param(y_dict,'background_distribution_branch','required')
                back_file = ROOT.TFile(back_file_name)
                back_tree = back_file.Get(back_tree_name)
                back_tree.Draw(back_branch_name+">>temp_hist","","goff")
                back_hist = ROOT.TH1F(ROOT.gDirectory.Get("temp_hist"))
                back_mag = back_hist.GetMean()

                renorm=1.0
                if(len(vars_dict)==2):
                    # The displayed fake data will be summed in the
                    # other dimension, so we need to rescale
                    other_var_num = (curr_var_num+1)%2
                    other_var = vars_dict[other_var_num]
                    other_n = read_param(other_var,'bins',30)
                    other_sig_shape = sum_weighted_shapes(signal_output_file,signals_dict,other_n,other_var_num,False,sig_weights)
                    other_back_shape = sum_weighted_shapes(back_output_file,back_dict,other_n,other_var_num,False,back_weights)
                    renorm = (sig_mag*np.sum(other_sig_shape)+back_mag*np.sum(other_back_shape))/(sig_mag+back_mag)
                
                y_arr = renorm*(sig_mag*tot_sig_shape+back_mag*tot_back_shape)
                if(read_param(y_dict,'print_mean_stddev',False)):
                    sig_err = sig_hist.GetStdDev()
                    fake_data_dict = read_param(prep_dict,'fake_data_settings','required')
                    text_to_print+="Sig=%.3f,Err=%.3f"%(sig_mag,sig_err)
                    true_sig_mag = read_param(fake_data_dict,'fake_signal_magnitude',0.0)
                    if(true_sig_mag!=0):
                        sig_pct = 100.0*sig_err/true_sig_mag
                        text_to_print+=",Pct=%.2f%%;"%sig_pct
                    else:
                        text_to_print+=";"
                    back_err = back_hist.GetStdDev()
                    text_to_print += "Back=%.3f,Err=%.3f"%(back_mag,back_err)
                    true_back_mag = read_param(fake_data_dict,'fake_background_magnitude',0.0)
                    if(true_back_mag!=0):
                        back_pct = 100.0*back_err/true_back_mag
                        text_to_print += ",Pct=%.2f%%;"%(back_pct)
                    else:
                        text_to_print+=";"
                if(err_bar_type=="extracted"):
                    # Propogate signal error and background error
                    sig_err = np.full(len(y_arr),sig_hist.GetStdDev())
                    back_err = np.full(len(y_arr),back_hist.GetStdDev())
                    def div0( a, b ):
                        # ignore divide by 0 error, set elt to 0
                        with np.errstate(divide='ignore', invalid='ignore'):
                            c = np.true_divide( a, b )
                            c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
                        return c
                    y_err = y_arr*np.sqrt(div0(sig_err,sig_mag*tot_sig_shape)**2+\
                                          div0(back_err,back_mag*tot_back_shape)**2 )
                sig_file.Close()
                back_file.Close()
            if(err_bar_type=="frac_sig"):
                y_err = err_bar_val*sig_mag*tot_sig_shape
            elif(err_bar_type=="frac_back"):
                y_err = err_bar_val*back_mag*tot_back_shape
            # For 2D data, the expected counts will have to be scaled according
            # to the normalization of the other dimension
            renorm = 1.0
        else:
            # Read y points from a root or R file
            file_name = read_param(y_dict,'file_name','required')
            file_type = read_param(y_dict,'file_type','required')
            tree_name = read_param(y_dict,'tree_name','required')
            if(file_type=="root"):
                branch_name = read_param(y_dict,'branch_name','required')
                y_arr = readTTree(file_name,tree_name,branch_name)
            elif(file_type=="R"):
                temp_R_dict = pystan.misc.read_rdump(file_name)
                y_arr = read_param(temp_R_dict,tree_name,'required')
                y_arr = np.array(y_arr)
                if(len(y_arr.shape)==2):
                    if(curr_var_num==1):
                        #y_arr = y_arr[0,:]
                        y_arr = np.sum(y_arr,0)
                    else:
                        #y_arr = y_arr[:,0]
                        y_arr = np.sum(y_arr,1)
                elif(len(y_arr.shape)>2):
                    print("ERROR: Fake data with more than 2 elts can't yet be plotted")
                    y_arr = np.zeros(n)
            else:
                print("%s is not a valid type, please select \"root\" or \"R\""%file_type)
        y_arr = np.array(y_arr)
        y_arr = y_arr.astype('float')
        # Calculate error bars if they were not calculated above
        if(err_bar_type=="frac_tot"):
            y_err = err_bar_val*y_arr
        elif(err_bar_type=="abs"):
            print("Absolute magnitude of error bars currently uses the given value for each bin, so the given value should be the bin width times the rate")
            y_err = np.full(n,err_bar_val)
        elif(err_bar_type=="poisson"):
            y_err = np.sqrt(y_arr)

        if(ymin==""):
            ymin = min(y_arr)
        else:
            ymin = min(ymin,min(y_arr))
        if(ymax==""):
            ymax = max(y_arr)
        else:
            ymax = max(ymax,max(y_arr))

        objs_to_print.append(ROOT.TGraphErrors(n,x_arr,y_arr,np.zeros(n),y_err))
        # End for loop
    ymin = ymin-0.1*(ymax-ymin)
    ymax = ymax+0.1*(ymax-ymin)
        
    txt_to_print_arr = []
    txt_xs = []
    txt_ys = []
    txt_xe = []
    txt_ye = []
    if(text_to_print!=""):
        txt_to_print_arr.append(text_to_print)
        txt_xs.append(0.01)
        txt_ys.append(0.01)
        xmax = min(len(text_to_print)*0.01+0.01,0.99)
        txt_xe.append(xmax)
        txt_ye.append(0.06)

    root_make_plot(objs_to_print,out_path,draw_opts,
                   title,label_x,label_y,
                   labels,lxs,lys,lxe,lye,leg_opts,
                   colors,[],marker_styles,
                   txt_to_print_arr,txt_xs,txt_ys,txt_xe,txt_ye,
                   ymin,ymax,set_fill_color=True)
    return
# ----------------------------------------------------------------------
def print_table(table_dict,out_dir,prep_dict=[]):
    root_files = read_param(table_dict,'files','required')
    trees = read_param(table_dict,'trees','required')
    branches = read_param(table_dict,'branches','required')
    print_means = read_param(table_dict,'print_means','required')
    print_errs = read_param(table_dict,'print_errs','required')
    print_pcts = read_param(table_dict,'print_pcts','required')
    
    
    out_files = []
    delimiters = []
    end_of_lines = []
    pct_signs = []
    
    if(read_param(table_dict,'do_tex_table',False)):
        tex_path = out_dir+read_param(table_dict,'tex_output_name','table.tex')
        temp_file = open(tex_path,"a+")
        out_files.append(temp_file)
        delimiters.append(" & ")
        end_of_lines.append("\\\\\n")
        pct_signs.append("\%")

    for i,f in enumerate(out_files):
        str = ""
        if(read_param(table_dict,'print_true_sig',False)):
            fake_data_dict = read_param(prep_dict,'fake_data_settings','required')
            true_sig_mag = read_param(fake_data_dict,'fake_signal_magnitude','required')
            str += "%.2f"%true_sig_mag+delimiters[i]
        if(read_param(table_dict,'print_true_back',False)):
            fake_data_dict = read_param(prep_dict,'fake_data_settings','required')
            true_back_mag = read_param(fake_data_dict,'fake_background_magnitude','required')
            str += "%.2f"%true_back_mag+delimiters[i]
        if(read_param(table_dict,'print_range',False)):
            var_num = read_param(table_dict,'var_range_num',0)
            vars_dict = read_param(prep_dict,'indep_vars','required')
            curr_var = vars_dict[var_num]
            lb = read_param(curr_var,'lower_bound','required')
            ub = read_param(curr_var,'upper_bound','required')
            r = ub-lb
            str += "%.1f"%r+delimiters[i]
        for j in range(0,len(root_files)):
            curr_file = ROOT.TFile(root_files[j])
            curr_tree = curr_file.Get(trees[j])
            curr_tree.Draw(branches[j]+">>temp_hist","","goff")
            curr_hist = ROOT.TH1F(ROOT.gDirectory.Get("temp_hist"))
            if(print_means[j]):
                str += "%.2f"%curr_hist.GetMean()+delimiters[i]
            if(print_errs[j]):
                str += "%.3f"%curr_hist.GetStdDev()+delimiters[i]
            if(print_pcts[j]):
                err = curr_hist.GetStdDev()
                if(j==0):
                    mag = true_sig_mag
                else:
                    mag = true_back_mag
                if(float(mag)!=0):
                    pct = 100.0*float(err)/float(mag)
                    str += "%.2f"%pct+pct_signs[i]+delimiters[i]
        str = str[:-len(delimiters[i])]
        str+=end_of_lines[i]
        f.write(str)
        f.close()
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
    plot_dict = read_param(cdata,'post_morpho_plots','required')
    if(read_param(plot_dict,'do_post_morpho_plots',True)):
        print("Starting plotting")
        # Set environment to make plots look nice
        sensitivitytools.root_tools.set_root_env()
        out_dir = read_param(plot_dict,'plots_output_directory','required')
        print("Plots will be stored in %s"%out_dir)
        create_path(out_dir,False)
        out_prefix = read_param(plot_dict,'plots_output_prefix','')
        title_postfix = read_param(plot_dict,'title_postfix','')
        if(read_param(plot_dict,'do_histograms',False)):
            print("Creating histogram plots")
            hist_dict = read_param(plot_dict,'histograms','required')
            for hist_settings in hist_dict:
                make_hist_plot(hist_settings,out_dir,out_prefix,title_postfix)
        if(read_param(plot_dict,'do_correlation_plots',False)):
            print("Creating correlation plots")
            corr_dict = read_param(plot_dict,'correlation_plots','required')
            for corr_settings in corr_dict:
                make_corr_plot(corr_settings,out_dir,out_prefix,title_postfix)
        if(read_param(plot_dict,'do_data_plots',False)):
            print("Creating data plots")
            dat_plot_dict = read_param(plot_dict,'data_plots','required')
            # if computed_points="true_data" or "extracted_data", the
            # preprocessing dictionary will be used to find shapes
            prep_dict = read_param(cdata,'preprocessing',[])
            for data_settings in dat_plot_dict:
                make_data_plot(data_settings,out_dir,out_prefix,
                               title_postfix,prep_dict)
        if(read_param(plot_dict,'do_print_table',False)):
            table_dict = read_param(plot_dict,'print_table','required')
            prep_dict = read_param(cdata,'preprocessing',[])
            print_table(table_dict,out_dir,prep_dict)
            
    print("Plotting finished")
