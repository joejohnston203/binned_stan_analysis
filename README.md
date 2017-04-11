# ricochet_sensitivity
Stan/morphomodel to analyze the sensitivity of Ricochet at the Double Chooz Reactor

Requirements
======

  Designed to work with version 1.1.5 of morpho. https://github.com/project8/morpho/tree/v1.1.5
  
  You will need to install root in order to have access to pyroot.
  
Directory Structure
======

  scripts

    The .yaml configuration files for stan models should be stored here. The
    .yaml file will also include parameters that will be used for preprocessing,
    such as generating fake data, under the header 'preprocessing:'.

    See rate_shape_analyzer.yaml for an example of how .yaml files should be
    written.

  ricfunctions

    Functions specific to ricochet, such as spectral shapes and time dependence
    shapes. The .yaml file should specifiy the files where the relevant
    functions are located.

    When adding signals or backgrounds to the model, you should only need to
    edit the .yaml config file and the functions in ricfunctions that the .yaml
    file references.
    
  preprocessing
  
    shape_fakedata_generator.py: Can generates the shape of signals and
    backgrounds using inputs from a .yaml file, then stores the shapes of
    the signals and backgrounds in the data directory. Informational
    outputs, such as plots of each background shape, can optionally be
    generated and stored. The generated shapes and inputs from a .yaml config
    file can then be used to generate fake data.

  data

    Information about signals and backgrounds, such as shapes and uncertanties,
    will be stored here. Generated fake data will also be stored here. Nothing
    in this folder should be edited by hand, because it should all be generated
    by scripts in the preprocessing folder.

  models

    Stan models (.stan files) are stored here. binned_data_analysis.stan is a
    very general stan model that analyzes binned data. The binned data can
    result from multiple signals and backgrounds, and can depend on multiple
    variables (such as recoil energy and time).

  results

    Outputs of the stan model will be stored here.

Running
======
