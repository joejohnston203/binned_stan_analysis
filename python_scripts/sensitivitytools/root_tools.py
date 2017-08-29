#!/usr/bin/env python
#
# Helper code to work with root files and make plots in root
# ----------------------------------------------------------------------
# Author: J. P. Johnston <jpj13@mit.edu>
#
# Some code is based on morpho
# morpho is located at https://github.com/project8/morpho/tree/v1.1.5
# The authors of morpho are:
#   J. A. Formaggio <josephf@mit.edu>
#   T. E. Weiss <tweiss@mit.edu>
#   M. Guigue <mathieu.guigue@pnnl.gov>
#   J. N. Kofron <jared.kofron@gmail.com>
# ----------------------------------------------------------------------

import ROOT as ROOT# import ROOT, TStyle, TCanvas, TH1F, TGraph, TLatex, TLegend, TFile, TTree, TGaxis, TRandom3, TNtuple, TTree
from ROOT import gStyle
from ROOT import gROOT

import numpy as np
from array import array
# ----------------------------------------------------------------------
# Set the style to make root plots look nice
def set_root_env():
  #//TStyle* genieStyle = new TStyle("genieStyle", "GENIE Style")
  #//set the background color to white
  gStyle.SetFillColor(10)
  gStyle.SetFrameFillColor(10)
  gStyle.SetCanvasColor(10)
  gStyle.SetPadColor(10)
  gStyle.SetTitleFillColor(0)
  gStyle.SetStatColor(10)
  
  #dont put a colored frame around the plots
  gStyle.SetFrameBorderMode(0)
  gStyle.SetCanvasBorderMode(0)
  gStyle.SetPadBorderMode(0)
  gStyle.SetLegendBorderSize(3)
  
  #use the primary color palette
  #gStyle.SetPalette(1,0)
  
  #set the default line color for a histogram to be black
  gStyle.SetHistLineColor(ROOT.kBlack)
  
  #set the default line color for a fit function to be red
  gStyle.SetFuncColor(ROOT.kRed)
  
  #make the axis labels black
  gStyle.SetLabelColor(ROOT.kBlack,"xyz")
  
  #set the default title color to be black
  gStyle.SetTitleColor(ROOT.kBlack)
  
  #set the margins
  gStyle.SetPadBottomMargin(0.18)
  gStyle.SetPadTopMargin(0.08)
  gStyle.SetPadRightMargin(0.08)
  gStyle.SetPadLeftMargin(0.17)
  
  #set axis label and title text sizes
  gStyle.SetLabelFont(42,"xyz")
  gStyle.SetLabelSize(0.04,"xyz")
  gStyle.SetLabelOffset(0.015,"xyz")
  gStyle.SetTitleFont(42,"xyz")
  gStyle.SetTitleSize(0.04,"xyz")
  gStyle.SetTitleOffset(1.4,"y")
  gStyle.SetTitleOffset(1.3,"x")
  gStyle.SetStatFont(42)
  gStyle.SetStatFontSize(0.07)
  gStyle.SetTitleBorderSize(1)
  gStyle.SetStatBorderSize(0)
  gStyle.SetTextFont(42)
  gStyle.SetTitleW(0.5)
  gStyle.SetTitleH(0.1)
  
  #set line widths
  gStyle.SetFrameLineWidth(2)
  gStyle.SetFuncWidth(2)
  gStyle.SetHistLineWidth(2)
  
  #set the number of divisions to show
  gStyle.SetNdivisions(506, "xy")
  #gStyle.SetPadTickX(-50202)
  
  #turn off xy grids
  gStyle.SetPadGridX(0)
  gStyle.SetPadGridY(0)
  
  #set the tick mark style
  gStyle.SetPadTickX(1)
  gStyle.SetPadTickY(1)
  
  #turn off stats
  gStyle.SetOptStat(0)
  gStyle.SetOptFit(0)
  
  #marker/line settings
  #gStyle.SetMarkerStyle(20)
  gStyle.SetMarkerSize(.95)#0.7
  gStyle.SetLineWidth(2) 
  gStyle.SetErrorX(0)
  gStyle.SetHistLineStyle(0) #It was 3 for a dotted line
  
  #done
  gStyle.cd()
  gROOT.ForceStyle()
#  gStyle.ls()
# ----------------------------------------------------------------------
# Add a title to a root plot
def add_plot_label(label,x,y,size=0.05,color=1,font=62,align=22):
  latex = ROOT.TLatex(x,y,label)
  latex.SetNDC()
  latex.SetTextSize(size)
  latex.SetTextColor(color)
  latex.SetTextFont(font)
  latex.SetTextAlign(align)
  latex.Draw()
  return latex
# ----------------------------------------------------------------------
# Creates a plot with all given objects and saves it
def root_make_plot(objs_to_plot,out_path,draw_opts=[],
                   title="",xlabel="",ylabel="",
                   legend_labels=[],lxs=0.7,lys=0.75,lxe=0.99,lye=0.9,
                   leg_opts=[],
                   colors=[],line_styles=[],marker_styles=[],
                   text_to_print=[],
                   txt_xs=[],txt_ys=[],txt_xe=[],txt_ye=[],
                   ymin="",ymax="",
                   titlex=0.5,titley=0.95,titlesize=0.035,
                   set_fill_color=False):
    c = ROOT.TCanvas()

    leg = ROOT.TLegend(lxs,lys,lxe,lye)
    for i in range(0,len(objs_to_plot)):
        if(i<len(colors)):
            objs_to_plot[i].SetLineColor(colors[i])
            objs_to_plot[i].SetMarkerColor(colors[i])
            if(set_fill_color):
              objs_to_plot[i].SetFillColor(colors[i])
        if(i<len(line_styles)):
            objs_to_plot[i].SetLineStyle(line_styles[i])
        if(i<len(marker_styles)):
            objs_to_plot[i].SetMarkerStyle(marker_styles[i])
        if(i<len(legend_labels) and legend_labels[i]!=""):
            if(i<len(leg_opts) and leg_opts[i]!=""):
                leg.AddEntry(objs_to_plot[i],legend_labels[i],leg_opts[i])
            else:
                leg.AddEntry(objs_to_plot[i],legend_labels[i])
            leg.Draw()
        if(i==0):
            objs_to_plot[i].GetXaxis().SetTitle(xlabel)
            objs_to_plot[i].GetYaxis().SetTitle(ylabel)
            if ymin!="" and ymax!="":
                objs_to_plot[i].GetYaxis().SetRangeUser(ymin,ymax)
            objs_to_plot[i].SetTitle("")
            if(i<len(draw_opts)):
                objs_to_plot[i].Draw(draw_opts[i])
            else:
                objs_to_plot[i].Draw()
            latex = add_plot_label(title,titlex,titley,size=titlesize)
        else:
            if(i<len(draw_opts)):
                objs_to_plot[i].Draw(draw_opts[i])
            else:
                objs_to_plot[i].Draw()

    # Add text
    pts = []
    for (i,txt) in enumerate(text_to_print):
        pts.append(ROOT.TPaveText(txt_xs[i],txt_ys[i],
                                  txt_xe[i],txt_ye[i],"NDC"))
        pts[i].AddText(txt)
        pts[i].Draw("same")

    c.SaveAs(out_path)
    
# ----------------------------------------------------------------------
# Create a histogram of the given branches
# file_names, tree_names, and branch_names all must be arrays of
#    the same length, giving the locations of the branches to plot
# labels: If len(labels)!=0, then a legend is added to the plot with
#         the given lables
# colors: 
def root_plot_histogram(file_names,tree_names,branch_names,
                        out_path,title="",x_label="",y_label="",
                        labels=[],leg_xstart=0.7,leg_ystart=0.75,
                        leg_xend=0.99,leg_yend=0.9,
                        colors=[],
                        quantiles_x_frac=[],print_mean_stddev=False):
    if not (len(file_names)==len(tree_names)
            and len(tree_names)==len(branch_names)):
        print("ERROR: file_name,tree_name, and branch_name arrays are different lengths")
        return

    # Make a canvas
    c = ROOT.TCanvas("c","c",800,600)
    
    if(print_mean_stddev):
        mean_err_str = ""

    file_objs = []
    tree_objs = []
    hist_objs = []
    draw_opts = []
    for i in range(0,len(file_names)):
        file_objs.append(ROOT.TFile(file_names[i],"curr_file"))
        tree_objs.append(file_objs[i].Get(tree_names[0]))
        tree_objs[i].Draw(branch_names[i]+">>temp_hist","","goff")
        hist_objs.append(ROOT.TH1F(ROOT.gDirectory.Get("temp_hist")))
        if(print_mean_stddev):
            mean = hist_objs[i].GetMean()
            err = hist_objs[i].GetStdDev()
            if(i<len(labels)):
               curr_label = labels[i]
            else:
               curr_label = branch_names[i]
            mean_err_str += curr_label + "=%.3f, err=%.3f; "%(mean,err)
        draw_opts.append("same")
    draw_opts[0]=""
        # End for loop
    txt_to_print = []
    txt_xs = []
    txt_ys = []
    txt_xe = []
    txt_ye = []
    if(print_mean_stddev):
        txt_to_print.append(mean_err_str[:-1])
        txt_xs.append(0.01)
        txt_ys.append(0.01)
        xmax = min(len(mean_err_str)*0.01+0.01,0.99)
        txt_xe.append(xmax)
        txt_ye.append(0.06)

    if(len(quantiles_x_frac)>0):
        q_str = "Quantiles: "
        for i in range(0,len(hist_objs)):
            if(len(hist_objs)>1):
                q_str+= "h%i- "%i
            q_x = np.asarray(quantiles_x_frac)
            q_y = np.zeros(len(q_x))
            hist_objs[i].GetQuantiles(len(q_x),q_y,q_x)
            for j in range(0,len(q_x)):
                q_str += "%.1f%%->%.3f, "%(100*q_x[j],q_y[j])
            q_str = q_str[:-2] + "; "
        txt_to_print.append(q_str[:-2])
        txt_xs.append(0.01)
        txt_ys.append(0.06)
        xmax = min(len(q_str)*0.008+0.01,0.99)
        txt_xe.append(xmax)
        txt_ye.append(0.1)

    root_make_plot(hist_objs,out_path,draw_opts,title,x_label,y_label,
                   labels,leg_xstart,leg_ystart,leg_xend,leg_yend,[],
                   colors,[],[],
                   txt_to_print,txt_xs,txt_ys,txt_xe,txt_ye)
    for f in file_objs:
        f.Close()
    return
# ----------------------------------------------------------------------
# Two branches, each of the same length, must be specified.
# A correlation plot of the first vs the second will then be saved.
def root_make_corr_plot(file_name_x,tree_name_x,branch_name_x,
                        file_name_y,tree_name_y,branch_name_y,
                        out_path,
                        label_x="",label_y="",
                        title=""):
    tc = ROOT.TChain()
    tc.SetName(tree_name_x)
    tc.Add(file_name_x)
    tc.SetName(tree_name_y)
    tc.Add(file_name_y)
    tc.Draw(branch_name_y+":"+branch_name_x+">>temp_hist","","goff")
    hist = ROOT.TH2F(ROOT.gDirectory.Get("temp_hist"))
    root_make_plot([hist],out_path,["","same"],title,label_x,label_y)
# ----------------------------------------------------------------------
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
