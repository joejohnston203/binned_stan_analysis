Stan/morpho models to analyze the sensitivity of any experiment where the data can be binned.

======

  Designed to work with version 1.4.1 of morpho. https://github.com/project8/morpho/tree/v1.4.1
  
  You will need to install root in order to have access to pyroot.

  The files from morpho_modules need to be copied into morpho/morpho before installing root.
  
Directory Structure
======

  config_builders

    Python scripts that generate a morpho config file and a stan model for a given analysis. For example, the cuore_config_builder.py script generates files necessary to do an analysis of 1-dimensional data. Using a config builder is NOT necessary, but it can be useful in an analysis that has a large number of parameters, such as the cuore analysis.

  docs

    Documentation. Needs to be updated.

  morpho_modules

    Morpho modules that need to be copied into morpho/morpho

  analyses

    Contains folders for multiple sample analyses. Inside each analysis folder, there generally are the following folders. Note however, that each of these names are defined in the yaml file used by morpho, and the names can be changed if they are also changed in the yaml file.

      scripts_config_builder

      .yaml configuration files that will be used by a config builder to make a
      Morpho script and a Stan model.

      scripts

        The .yaml configuration files for morpho should be stored here. The
        .yaml file will include parameters that will be used by morpho for preprocessing,
	running stan, and making plots.

      models

        Stan models (.stan files) are stored here.

      functions

        Functions used by the stan model should be stored here.

      data

        Any other inputs should be stored here, such as root files with the
	spectral shape of a parameter.

      The following folders should initially be empty when an analysis starts running,
      and they will then be populated by morpho:

      morpho_data

        Information about signals and backgrounds, such as shapes and uncertanties,
        will be stored here. Generated fake data will be stored here. Diagnostic outputs
        and plots, such as plots of signal and background shapes, can also be stored.

      results

        Outputs will be stored here. Specifically, stan will store a .root file with
        distributions of parameters, a .pkl file, and a cache_name_file.

      cache

        cached stan models will be stored here, so the model does not have to be recompiled
        every time it is run.

CUORE 0 Analysis
======

All analysis scripts should be run from the folder analyses/cuore0_analysis.

The cuore-0 analysis creates a model that is the sum of many parameters, with each paraeter representing some background in the cuore detector. It allows for multiple data sets, such as an M1 (single detector deposit) and an M2 (coincident deposits in two detectors) spectrum. For each parameter+data set combination, a spectrum for the parameter must be provided.

The cuore config builder can be used to create the yaml file. cuore0_config_builder_simplified.yaml is a simplified file that referes only to shapes that are stored in the github repository. It can be used to generate a morpho config file and a stan model with the command:

python ../../config_builders/cuore_config_builder.py -c scripts/config_builder/cuore0_config_builder_simplified.yaml

This will create a morpho scripts, "scripts/cuore0_analysis_simplified.yaml", and a stan model, "models/cuore0_analysis_simplified.stan".

This simplified model can then be run with morpho:

morpho -c scripts/cuore0_analysis_simplified.yaml

This will do the following:
  - Obtain the spectral shape for each parameter and store the shapes to an R file
  - Either generate fake data or access stored data, and store it to an R file
  - Run a stan model that fits for the magnitude of each parameter
  - Generate plots of the spectra and other useful quantities


Ricochet Analysis
======

The ricochet_analysis folder contains scripts scripts that were used to predict
the sensitivity of the Ricochet experiment for the case of 1D (time domain only)
and 2D (time+energy) data. It currently is NOT set up to use the morpho_modules,
and needs to be updated. Also, scripts currently are not working, and will need
some debugging to be fixed. The way it currently is intended to work is as follows:

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
