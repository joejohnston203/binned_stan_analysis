Stan/morpho model to analyze the sensitivity of any experiment where the data can be binned. Specific scripts analyze the sensitivity of Ricochet at the Double Chooz Reactor.

The frequentist_analysis folder contains a python program to perform a frequentist analysis. It takes configuration scripts for the Bayesian code and performs a frequentist analysis. See the end of this README for more details on its structure and use.

Requirements
======

  Designed to work with version 1.1.5 of morpho. https://github.com/project8/morpho/tree/v1.1.5
  
  You will need to install root in order to have access to pyroot.

  If you want to use the given example analysis, update activate_venv.sh to point to your virtualenv directory where morpho was installed.
  
Directory Structure
======

  The following directories contain files that will be copied and altered when creating new models:

  scripts

    The .yaml configuration files for stan models should be stored here. The
    .yaml file will include parameters that will be used for preprocessing,
    such as generating fake data, under the header 'preprocessing:'. It will also
    include parameters for plots, under the header 'post_morpho_plots:'.

    See example_ricochet_rate_analyzer.yaml or example_ricochet_rate_spec_analyzer.yaml
    for examples of how .yaml files should be created.

  ric_functions

    Functions and .csv files specific to the Ricochet experiment, such as spectral
    shapes and time dependence shapes. The .yaml file should specifiy the paths
    where the relevant files are located.

  models

    Stan models (.stan files) are stored here. For example,
    binned_data_analysis_1D_1back.stan is a model that takes the inputs from
    pre_morpho_processing.py, then runs an analysis of the generated 1-dimensional
    data assuming that it results from a single signal and a single background.
    binned_data_analysis_2D_4back_tot_bin_var.stan is a model that analyzes
    2-dimensional data assuming that it results from one signal and four backgrounds.
    This model also allows the total number of counts in each bin to vary by some amount
    (if such variation was turned on when preprocessing was done). (The tot_bin_var
    version of this model should only be used if the model actually does allow the
    counts to fluctuate in each bin because it does run more slowly).

    When creating a new analysis, a new model should be created that properly fits
    for the correct number of signals and backgrounds, as well as any gaussian
    fluctuations in the signals, backgrounds, and counts. The stan models should be
    fairly general. For example, binned_data_analysis_1D_2back.stan can be used
    to analyze any 1-dimensional data that results from 1 signal and 2 backgrounds,
    and allows the signal to gaussian fluctuate by some fraction in each bin, and/or
    to fluctuate globally.

  bash_scripts

    Scripts to iterate over configurations and run the analysis for each. For example,
    example_run_1D_analysis.sh copies the example_ricochet_rate_analyzer.yaml script,
    then edits the copy to run the analysis with different signal and background
    magnitudes.

  The following directory should not be altered when creating new models:

  python_scripts
  
    pre_morpho_processing.py: Can generates the shape of signals and
    backgrounds using inputs from a .yaml file, then stores the shapes of
    the signals and backgrounds in the data directory. Informational
    outputs, such as plots of each background shape, can optionally be
    generated and stored. The generated shapes and inputs from a .yaml config
    file can then be used to generate fake data.

    post_morpho_plots.py: Creates histograms of analysis parameters, correlation
    plots, and plots of spectra or rate vs time. Takes inputs from a .yaml file.

    root_tools.py: Contains tools to read and write root files, as well as make
    plots using root.
    dict_tools.py: Contains tools to read a dictionary from a .yaml file, and
    work with that dictionary.

    run_analysis.sh: Script to run preprocessing, stan, and make plots for
    a given .yaml file.

  The sample analysis will create the following three directories while running.
  Sample outputs from ./bash_scripts/example_simple_run.sh are currently stored
  here, but these directories could be deleted before running an analysis. The
  locations of these directories can be edited in the .yaml file:

  data

    Information about signals and backgrounds, such as shapes and uncertanties,
    will be stored here. Generated fake data will be stored here. Diagnostic outputs
    and plots, such as plots of signal and background shapes, can also be stored. Nothing
    in this folder should be edited by hand, because it will all be generated
    by scripts in the python_scripts folder.

  results

    Outputs will be stored here. Specifically, stan will store a .root file with
    distributions of parameters, a .pkl file, and a cache_name_file. post_morpho_plots.py
    will then store plots and a table of the results.

  cache

    cached stan models will be stored here, so the model does not have to be recompiled
    every time it is run.

Running
======

  Once you have created a yaml_script, the following commands run the analysis and create plots:

  python python_scripts/pre_morpho_processing.py -c $yaml_script
  
  morpho -c $yaml_script
  
  python python_scripts/post_morpho_plots.py -c $yaml_script

  To run the rate only ricochet analysis, use the script ./bash_scripts/example_simple_run.sh.
  Edit the variable virtualenv_script to point to the activation script for the
  virtualenv environment that you set up when installing stan. Then run the script from
  the top directory. This will use ./scripts/example_ricochet_rate_analyzer.yaml to
  do preprocessing, run stan, and create plots with post_morpho_plots.py

  example_1D_analysis.sh and example_2D_analysis.sh both iterate over different signal and
  background magnitudes, doing preprocessing, running stan, and creating plots for each
  signal and background combination. In order to get these to run, you should just need to
  update the variable virtualenv_script. However, these do have other variables that
  can be altered to change the runs. Note that these copy the .yaml file before editing it,
  so these can safely be run without altering the initial .yaml file.

  ./bash_scripts/ricochet_analysis/run_full_analysis.sh will sequentially run all
  configurations that went into the paper "Coherent Neutrino Scattering with Low 
  Temperature Bolometers at Chooz Reactor Complex." Note that each configuration is
  run many times with a single chain. This is so you can check the output plots to
  discard any chains that did not properly converge, then average the results of the
  remaining chains to get final results.

  The results that went into the Ricochet at Double Chooz paper can be accessed at:
  https://www.dropbox.com/sh/owa9y6qtmdcwenk/AABiqeQbo0WBhvnxuDFqY1VEa?dl=0

Frequentist Analysis
======

This directory contains code to do a frequentist analysis. It takes a .yaml configuration file, and uses the same pre_processing and post_morpho_plots sections as the .yaml files for a bayesian analysis. Instead of the morpho and stan sections, a section labelled frequentist contains settings for the frequentist analysis.

# Directory Structure (inside frequentist_analysis)

The following directories can be changed when creating new models

models

  Likelihood functions are stored here.

../scripts, ../ric_functions

  .yaml scripts for the frequentist analysis will be formatted very similarly
  to the scripts for the bayesian analysis. In fact, the preprocessing section
  can be identical. However, a frequentist section must be added with parameters
  in place of the stan and morpho sections, and the plotting section will be
  different.

The following should not be edited

analysis_scripts

  This contains the actual python code that does the analysis

data, results

  These folders are the currently specified locations in the .yaml script. If
  the analysis is called from inside frequentist_analysis, then these
  directories will also be created inside frequentist_analysis