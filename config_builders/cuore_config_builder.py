#======================================================
# generate_scripts.py
#
# Author: J. Johnston
# Date: Mar. 1, 2018
#
# Generate the yaml script and Stan model required by
# Morpho
#=======================================================

"""Generate a morpho config file and stan model for binned analysis

Many parameters, such as the locations where files are stored, are
shared between morpho modules. This module is designed such that
each piece of info is only input once, and is then carried through
to all morpho modules.

Classes:
  - StringBuilder: Class to generate a large string with indents
  - BinnedConfigBuilder: Class to generate morpho config and stan model

Functions:
  - generate_binned_config: Generate the config and stan model
  - __main__: Call generate_binned_config

ToDo:

I need to be able to multiply multiple likelihoods, eg M1+M2
spectrum. But I think that can be handled without changing
the processors I have written. I will still need a shape for
every parameter. I can call the fake data generator multiple
times to make multiple data sets. I can call the plotting
module multiple times, once for each data set.

(Note- I think it is safe to assume that all parameters are
rates or exposures, so the likelihood is the product of the
likelihood for each dataset, where for a given data set,
expected counts is the sum of the parameter magnitude
times the parameter shape)

Determine if there's a better way to pass the settings for
reading a file. Having fields for 'columns', 'var_names',
etc is a lot, it would be nice to just have one field.

Figure out functions folder. Currently I have it so you can
specify different functions files for each shape, which
is not right.
"""

import os
from collections import OrderedDict

import logging
logger = logging.getLogger()

import yaml
yaml.Dumper.ignore_aliases = lambda *args : True
from argparse import ArgumentParser

import numpy as np

from morpho.utilities.reader import read_param
#from morpho.utilities.file_reader import *
#from morpho.utilities.file_writer import *

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

class BinnedConfigBuilder:
    """Generate a Morpho config file and Stan model

    Use Config(params) to read in all inputs from a dictionary,
    then call Run() to generate the Morpho config and Stan model.
    params is a dictionary specifying all configuration for the
    analysis. It will contain the following keys:

    Args:
        generate_morpho_config: Boolean specifying whether the Morpho
            config file should be generated (Default=True)
        morpho_config_output_path: Output path for the Morpho config
            file (Default="scripts/binned_analysis.yaml")
        generate_stan_model: Boolean specifying wether the Stan model
            should be generated (Default=True)
        stan_model_output_path: Output path for the Stan model
            (Default="models/binned_analysis.stan")

        output_paths: Dictionary containing output paths and formats for
            the analysis. Should contain the following keys:
          - shape_output_dir: Output directory used to store parameter
            shapes (Default="data/binned_analysis/shapes")
          - misc_config_output_dir: Output directory used to store
            miscellaneous configuration for the Stan model, such as
            the number of bins.
            (Default="data/binned_analysis/misc_config")
          - fake_data_output_dir: Output directory used to store fake data
            if generate_fake_data is True.
            (Default="data/binned_analysis/")
          - morpho_output_dir: Output directory used to store the
            outputs of morpho (Default="results")
          - morpho_output_file: Filename used to store the
            outputs of morpho (Default="analysis_results.root")
          - morpho_output_tree: Tree used to store the results
            to file. (Default="analysis_parameters")
          - plots_output_dir: Output path used to store plots
            (Default="results/plots/")
          - plots_output_format: File extension of the output format
            (Default="png")
          - debug_dir: Directory used to store debugging output
            (Default="data/binned_analysis/debug")

        stan: Dictionary containing configuration for Stan
            (required). All elements will be placed in the Morpho
            config file. Should contain the following fields:
          - function_files_location: String specifying directory
            containing funtions (Default=None)
          - cache: String specifying location to store the cache
            (Default="./cache")
          - run: Dictionary containing configuration for the Stan
            sampler. Should contain the following keys:
              - algorithm: Sampling algorithm (Default="NUTS")
              - iter: Number of iterations (required)
              - warmup: Number of iterations to discard as warmup
                (required)
              - chain: Number of chains to run (Default=1)

        dimension: Dictionary defining the x dimension.
            (required). Must contain the following keys:
          - name: Name of the dimension (Will be used to name files and
            variables, so should not start with a number or symbol, and
            be careful of clashes with Stan functions, etc).
            (required).

        generate_fake_data: Boolean specifying whether fake data should
            be generated. If true, fake data will be generated, stored,
            then passed to all future modules. If False, then data will
            be loaded from file and passed to all future modules.
            (Default=False)
        num_data_sets: Number of data sets. (Default=1)
        data: List of dictionaries, with one for each data set. Each
            parameter will be assigned to one data set. Each data set
            should contain the following fields:
          - name: Name of the data set, used for storing and titling plots
            (Default="Data_Set_i")
          - load_data_path: Path to data that will be loaded if
            generate_fake_data is False.
            (required if generate_fake_data==False)
          - load_data_format: Format for loaded data. Options are "text",
            "root", or "R". (required if generate_fake_data==False)
          - load_data_variable: Variale used to access the data. For root,
            should be a length 2 list specifying a tree and a branch name.
            For R, should be a string specifying the variable name.
            (required if generate_fake_data==False)
          - binning: Settings for determining binning. The following
            should be included for each data set (required)
              - min, max: Minimum and maximum. Optional if a rebin region
                is given, but if they are given, they will override all
              - binning described below if it would go outside the range.
              - min_bin_width: Bin width to be used when creating binning
                between max and min. The largest number of bins possible will
                be created such that the bin width is still larger than the minimum.
              - include_peaks: Whether bins should be added for the given peaks
              - peak_means_path, peak_means_format, peak_means_variables,
                peak_widths_path, peak_widths_format, peak_widths_variables:
                Files with lists of peaks to be included in the binning, such
                that the binning still respects the given minimum bin width
              - merge_low_stats_bins: If data is given, then the number of counts
                in each bin can be considered, and bins can be merged if the
                number of counts is less than a given min
              - min_bin_counts: Minimum number of counts if merge_low_stats_bins
                is true
              - rebin_regions: Dictionaries giving a path, format, and variables
                with lists of bin edges for regions that should be overridden

        parameters: List of dictionaries containing the parameters that
            will go into the model. Each should contain the following
            keys:
          - name: Name of the parameter (required). 
            (Will be used to name files and
            variables, so should not start with a number or symbol, and
            be careful of clashes with Stan functions, etc).
          - lower_bound: Lower bound in Stan (Default=0.)
          - upper_bound: Upper bound in Stan (Required)
          - prior: String, in Stan code, with the prior. Will add
            'p_name ~ prior_string' to the morpho model. (Default=None)
          - hierarchical_gaussian: String specifying whether this
            parameter should be treated as a hierarchical parameter with
            the parameter being allowed to fluctuate about the global
            value defined by the prior. Options are "False", "bin",
            or "dimension 0". "False" does no fluctuation. "bin" allows
            the parameter to fluctuate independently in every bin.
            "dimension 0" allows the parameter to fluctuate in each bin
            of the 0th dimension (where any dimension can be specified).
            (Default="False")
          - hierarchical_guassian_fraction: Fraction specifying the
            fractional fluctuation (Default=0.)
          - hierarchical_gaussian_center: Whether the hierarchical
            parameters for this parameter should be centered at 0.
          - fake_data_magnitude: Magnitude used to generate fake data.
            (required if generate_fake_data is True)
          - shapes: List of dictionaries, one dictionary
            for each data set defined above.
            Each dictionary should contain the following keys:
              - path: Path to the file containing the shape
              - renormalize: Boolean specifying whether the shape should
                be renormalized to 1 (Default=False)
              - multiply_shape: Float that multiplies the shape. If
                renormalize is True, then the shape is multiplied after
                it is renormalized to 1. (Default=1.0)
              - number_save_type: String specifying way to save number,
                from numpy data types (eg 'int64' or 'float64').
                (Default='float64').
              - format: Format of the file containing the shape.
                options are "text values", "text counts", "text function",
                "root values", "root counts", "root function", "R values",
                "R counts", "R function", or "python function".
                "text", "root", "R", or "python" specifies the filetype.
                "values" specifies that there will be a 1D list of values to
                be histogrammed. "counts" specifies that the number of counts
                in each bin will be specified (assuming bins from binning_file
                or n_bins). "function" specifies that there will be two
                columns of points defining the spectrum. In the case of
                function, the integral over the bin will be stored as the shape.
              - columns: List of column(s) including the relevant data. One column
                number should be included for "values" or "counts", two for "function".
                Required  if format includes "text" or "numpy". 
              - tree: Root tree name used to access data. Required if format
                includes "root".
              - branches: list of root branch name(s) used to access data. 
                (See columns). Required if format includes "root".
              - var_names: Variable names in the "R" file. (See columns).
                Required if format includes "R".
              - module: Module where the python function will be found.
                Required if format includes "python". In this case, path
                should specify the path to the directory where the module
                is located.
              - function: Function name. Required if format includes "python".
              - function_options: Dictionary with any named arguments to pass to
                the function. (Default={})

        plot: Dictionary containing configuration for the output plots.
            Should contain the following keys:
          - param_histos: Dictionary containing configuration for plots
            containing a histogram for each parameter.
              - make_plots: Boolean specifying whether parameter
                histogram plots should be made (Default=True)
              - output_dir: Output directory inside the plots folder
                (Default="param_histos")
          - correlations: Dictionary containing configuration for
            plots of the correlations between parameters
              - plot_aposteriori_distributions: Boolean specifying
                whether a grid of all parameter combinations should be
                made with a 2D histogram of each combination.
              - plot_correlation_factors: Boolean specifying whether
                a grid of all parameter combinations should be made
                with color representing the correlation factor between
                each combination.
          - binned_spectra: Dictionary containing configuration for
            plots of the spectra. If the analysis were multidimensional,
            then each dimension would be plotted while limiting to the
            first bin of the other dimensions.
              - make_plots: Boolean specifying whether plots of
                the binned spectra should be made
              - output_dir: Output directory inside the plots folder
                (Default="spectra")
          - which_plot: List of dictionaries containing any additional
            plots that should be made. Will be appended to the which_plot
            section of the Morpho config file. (Default=[])

    Returns:
        None: The morpho config file and the stan model are output.
    """
    def __init__(self, name, *args, **kwargs):
        self.__name = name
        return

    def Configure(self, params):
        self.generate_morpho_config = \
            read_param(params, 'generate_morpho_config', True)
        self.morpho_config_output_path = \
            read_param(params, 'morpho_config_output_path',
                       "scripts/binned_analysis.yaml")
        self.generate_stan_model = \
            read_param(params, 'generate_stan_model', True)
        self.stan_model_output_path = \
            read_param(params, 'stan_model_output_path',
                       "models/binned_analysis.stan")

        paths_dict = read_param(params, 'output_paths', {})
        self.shape_output_dir = \
            read_param(paths_dict, 'shape_output_dir',
                       "data/binned_analysis/shapes")
        self.misc_config_output_dir = \
            read_param(paths_dict, 'misc_config_output_dir',
                       "data/binned_analysis/misc_config")
        self.fake_data_output_dir = \
            read_param(paths_dict, 'fake_data_output_dir',
                       "data/binned_analysis/")
        self.morpho_output_dir = \
            read_param(paths_dict, 'morpho_output_dir', "results")
        self.morpho_output_file = \
            read_param(paths_dict, 'morpho_output_file', "analysis_results.root")
        self.plots_output_dir = \
            read_param(paths_dict, 'plots_output_dir', "results/plots")
        self.plots_output_format = \
            read_param(paths_dict, 'plots_output_format', "png")
        self.debug_dir = \
            read_param(paths_dict, 'debug_dir',
                       "data/binned_analysis/debug")

        self.morpho_output_tree = \
            read_param(paths_dict, 'morpho_output_tree', "analysis_parameters")

        self.stan_dict = read_param(params, 'stan', {})
        self.stan_dict["model"] = read_param(self.stan_dict, "model", "required")
        self.stan_dict["model"]["cache"] = \
            read_param(self.stan_dict, 'model.cache', "./cache")

        if self.generate_stan_model:
            self.stan_dict["model"]["file"] = \
                self.stan_model_output_path
        else:
            self.stan_dict["model"]["file"] = \
                read_param(self.stan_dict["model"], "file", "required")
        self.stan_dict["run"] = read_param(self.stan_dict, 'run', {})
        self.stan_dict["run"]["algorithm"] = \
            read_param(self.stan_dict["run"], 'algorithm', "NUTS")
        self.stan_dict["run"]["iter"] = \
            read_param(self.stan_dict["run"], 'iter', 'required')
        self.stan_dict["run"]["warmup"] = \
            read_param(self.stan_dict["run"], 'warmup', 'required')
        self.stan_dict["run"]["chain"] = \
            read_param(self.stan_dict["run"], 'chain', 1)

        self.dim_name = read_param(params, 'dimension_name', 'required')

        self.generate_fake_data = read_param(params, 'generate_fake_data', False)
        self.data_sets = read_param(params, 'data', 'required')
        self.num_data_sets = read_param(params, 'num_data_sets', 1)
        self.data_set_names = []
        self.load_data_paths = []
        self.load_data_formats = []
        self.load_data_variables = []
        self.binnings = []
        for i_data,data in enumerate(self.data_sets):
            self.data_set_names.append(read_param(data, 'name', "Data_Set_%i"%i_data))
            self.load_data_paths.append(read_param(data, 'load_data_path', None))
            self.load_data_formats.append(read_param(data, 'load_data_format', None))
            if not self.generate_fake_data and not self.load_data_formats[-1]=="root values":
                logger.error("Currently 'root values' is the only supported format for data."
                             + "%s' is invalid."%self.load_data_formats[-1])
            self.load_data_variables.append(read_param(data, 'load_data_variable', None))
            self.binnings.append(read_param(data, 'binning', 'required'))

        self.shapes_files = []
        self.binning_files = []
        for i_data in range(self.num_data_sets):
            self.shapes_files.append(self.shape_output_dir + "/" +
                                     self.data_set_names[i_data] + "_shapes.out")
            self.binning_files.append(self.shape_output_dir + "/" +
                                      self.data_set_names[i_data] + "_binnings.out")
        self.fake_data_file = self.fake_data_output_dir + "/fake_data.out"

        self.parameters = \
            read_param(params, 'parameters', 'required')

        self.num_params = len(self.parameters)
        self.param_names = []
        self.param_lower_bounds = []
        self.param_upper_bounds = []
        self.param_priors = []
        self.param_hierarchical_gaussians = []
        self.param_hierarchical_gaussian_fractions = []
        self.param_hierarchical_gaussian_centers = []
        self.fake_data_magnitudes = []
        self.param_shapes = [] 

        for p in self.parameters:
            self.param_names.append(
                read_param(p, 'name', 'required'))
            self.param_lower_bounds.append(
                read_param(p, 'lower_bound', 0.))
            self.param_upper_bounds.append(
                read_param(p, 'upper_bound', 'required'))
            self.param_priors.append(
                read_param(p, 'prior', None))
            self.param_hierarchical_gaussians.append(
                read_param(p, 'hierarchical_gaussian', False))
            if (isinstance(self.param_hierarchical_gaussians[-1], basestring) and
                "dimension" in self.param_hierarchical_gaussians[-1]):
                self.param_hierarchical_gaussians[-1] = \
                    self.param_hierarchical_gaussians[-1].split()[1]
            else:
                self.param_hierarchical_gaussians[-1]= \
                    [self.param_hierarchical_gaussians[-1]]
                
            self.param_hierarchical_gaussian_fractions.append(
                read_param(p, 'hierarchical_gaussian_fraction', 0.0))
            self.param_hierarchical_gaussian_centers.append(
                read_param(p, 'hierarchical_gaussian_center', 0.0))
            if self.generate_fake_data:
                self.fake_data_magnitudes.append(
                    read_param(p, 'magnitude', 'required'))

            self.param_shapes.append(read_param(p, 'shapes', 'required'))
            for s in self.param_shapes[-1]:
                s["renormalize"] = read_param(s, 'renormalize', False)
                s["multiply_shape"] = read_param(s, 'multiply_shape', 1.0)

        plot_dict = read_param(params, 'plot', {})
        param_histo_dict = read_param(plot_dict, 'param_histos', {})
        self.plot_param_histos = \
            read_param(param_histo_dict, 'make_plots', True)
        self.output_dir_param_histos = \
            read_param(param_histo_dict, 'output_dir', "param_histos")
        correlations_dict = read_param(plot_dict, 'correlations', {})
        self.plot_aposteriori_distributions = \
            read_param(correlations_dict,
                       'plot_aposteriori_distributions', True)
        self.plot_correlation_factors = \
            read_param(correlations_dict,
                       'plot_correlation_factors', True)
        binned_spectra_dict = read_param(plot_dict, 'binned_spectra', {})
        self.plot_binned_spectra = \
            read_param(binned_spectra_dict, 'make_plots', True)
        self.output_dir_binned_spectra = \
            read_param(binned_spectra_dict, 'output_dir', "spectra")
        self.which_plot = read_param(plot_dict, 'which_plot', [])
        return

    def get_morpho_config(self):
        """Build a dictionary for morpho configuration

        Args:
            None: Config() must have been run

        Returns:
            string: Returns yaml dictionary that can be used to
            configure Morpho
        """
        # PyYAML always sorts dictionaries before dumping, but we do not
        # want to sort, because the order of the preprocessing matters
        class UnsortableList(list):
            def sort(self, *args, **kwargs):
                pass
        class UnsortableOrderedDict(OrderedDict):
            def items(self, *args, **kwargs):
                return UnsortableList(OrderedDict.items(self, *args,
                                                        **kwargs))
        yaml.add_representer(UnsortableOrderedDict,
                             yaml.representer.SafeRepresenter.represent_dict)
        
        morpho_config_dict = UnsortableOrderedDict()

        # Morpho global configuration
        morpho = UnsortableOrderedDict()
        morpho["do_preprocessing"] = True
        morpho["do_stan"] = True
        morpho["do_plots"] = True

        morpho_config_dict["morpho"] = morpho

        # Preprocessing configuration
        preprocessing = []

        bin_shapes_dicts = []
        for i_data in range(self.num_data_sets):
            bin_shapes = UnsortableOrderedDict()
            bin_shapes["module_name"] = "binned_shapes"
            bin_shapes["method_name"] = "generate_shapes"
            bin_shapes["generate_shapes"] = True
            bin_shapes["output_dir"] = self.shape_output_dir
            bin_shapes["output_path_prefix"] = self.data_set_names[i_data] + "_"
            bin_shapes["output_format"] = "R"
            bin_shapes_dimension = UnsortableOrderedDict()
            bin_shapes_dimension["name"] = self.dim_name
            bin_shapes_dimension["binning"] = self.binnings[i_data]
            bin_shapes["dimensions"] = [bin_shapes_dimension]
            bin_shapes_params = list()
            for i_param in range(self.num_params):
                bin_shapes_params.append(
                    {
                        "name":"%s_%i"%(self.param_names[i_param],i_data),
                        "regenerate":True,
                        "shapes":[self.param_shapes[i_param][i_data]]
                    }
                )
            # Put the shape of the data in the file
            if not self.generate_fake_data:
                bin_shapes["binning_data_path"] = self.load_data_paths[i_data]
                bin_shapes["binning_data_format"] = self.load_data_formats[i_data]
                bin_shapes["binning_data_variables"] = self.load_data_variables[i_data]
                if self.load_data_formats[i_data]=="root values":
                    if len(self.load_data_variables[i_data])<3:
                        self.load_data_variables.append("") # Third element should be a cut
                    bin_shapes_params.append(
                        {
                            "name":self.data_set_names[i_data],
                            "regenerate":True,
                            "shapes": [{
                                "path": self.load_data_paths[i_data],
                                "format": self.load_data_formats[i_data],
                                "tree": self.load_data_variables[i_data][0],
                                "branches": [self.load_data_variables[i_data][1]],
                                "cut": self.load_data_variables[i_data][2],
                                "renormalize":False,
                                "multiply_shape":1.0,
                                "number_save_type": "int64"
                            }]
                        }
                    )
                else:
                    logger.error("Currently 'root values' is the only supported data format")
            bin_shapes["parameters"] = bin_shapes_params
            bin_shapes_dicts.append(bin_shapes)

        fake_data_dicts = []
        if self.generate_fake_data:
            for i_data in range(self.num_data_sets):
                fd_dict = UnsortableOrderedDict()
                fd_dict["module_name"] = "binned_fake_data"
                fd_dict["method_name"] = "binned_fake_data"
                fd_dict["output_dir"] = self.fake_data_output_dir
                fd_dict["output_path_prefix"] = "fake_data"
                fd_dict["output_format"] = "R"
                fd_dict["output_variable"] = "%s_%s"%(self.data_set_names[i_data],
                                                        self.dim_name)
                fd_params = list()
                for i_param, p in enumerate(self.parameters):
                    fd_params.append({})
                    fd_params[-1]["name"] = self.param_names[i_param]
                    fd_params[-1]["magnitude"] = self.fake_data_magnitudes[i_param]
                    fd_params[-1]["shapes"] = [{
                        "path": self.shapes_files[i_data],
                        "format":"R",
                        "variables":"%s_%i_%s"%(self.param_names[i_param],
                                                i_data, self.dim_name)
                    }]
                fd_dict["parameters"] = fd_params
                fake_data_dicts.append(fd_dict)

        which_prep = UnsortableOrderedDict()
        which_prep["which_pp"] = bin_shapes_dicts+fake_data_dicts
        morpho_config_dict["preprocessing"] = which_prep

        # Stan configuration
        stan_data_files = []
        for i_data in range(self.num_data_sets):
            stan_data_files.append(
                {"name": self.shapes_files[i_data],
                 "format":"R"})
            stan_data_files.append(
                {"name":self.binning_files[i_data],
                    "format":"R"})
        if self.generate_fake_data:
            stan_data_files.append(
                {"name":self.fake_data_file,
                 "format":"R"})

        self.stan_dict["data"] = {
            "files": stan_data_files,
            "parameters":[
                {
                    
                }
            ]
        }
        params_to_save = []
        for i_param in range(len(self.parameters)):
            params_to_save.append(
                {"variable": "rate_%s"%self.param_names[i_param],
                 "root_alias": "rate_%s"%self.param_names[i_param]}
            )
        self.stan_dict["output"] = {
            "name":self.morpho_output_dir+"/"+self.morpho_output_file,
            "format":"root",
            "tree":self.morpho_output_tree,
            "save_cache_name":self.misc_config_output_dir+"/cache_name_file.txt",
            "fit": self.morpho_output_dir+"/analysis_fit.pkl",
            "branches":params_to_save
        }
        morpho_config_dict["stan"] = self.stan_dict

        # Plotting configuration
        plotting = []

        binned_spectra_dict = UnsortableOrderedDict()
        binned_spectra_dict["module_name"] = "binned_spectra"
        binned_spectra_dict["method_name"] = "reconstructed_spectrum"
        binned_spectra_dict["output_dir"] = self.plots_output_dir
        binned_spectra_dict["individual_param_output_dir"] = self.plots_output_dir +\
                                                             "/param_reconstructions"
        binned_spectra_dict["output_path_prefix"] =  ""
        binned_spectra_dict["output_format"] = self.plots_output_format
        binned_spectra_dict["store_param_dists"] = True
        binned_spectra_dict["store_param_fractions"] = False
        binned_spectra_dict["store_param_dists_dir"] = self.plots_output_dir+"/param_dists"
        binned_spectra_dict["make_individual_spectra"] = False
        binned_spectra_dict["make_stacked_spectra"] = False
        binned_spectra_dict["make_unstacked_spectra"] = False
        binned_spectra_dict["make_reconstruction_plot"] = False
        binned_spectra_dict["make_diff_plot"] = False
        binned_spectra_dict["make_residual_pull_plot"] = False
        binned_spectra_dict["make_data_model_ratio_plot"] = False
        binned_spectra_dict["make_chi2_vs_dof_plot"] = False
        binned_spectra_dict["make_data_plot"] = False
        binned_spectra_dict["make_param_dist_plots"] = True
        binned_spectra_dict["binning_file"] = self.binning_files[0]
        binned_spectra_dict["binning_file_format"] = "R"
        binned_spectra_dict["binning_file_variable"] = self.dim_name
        binned_spectra_dict["divide_by_bin_width"] = True
        binned_spectra_dict["xlabel"] = "Energy (keV)"
        binned_spectra_dict["ylabel"] = "Counts/keV"
        binned_spectra_dict["ylog"] = True
        binned_spectra_dict["title_prefix"] = ""
        if self.generate_fake_data:
            binned_spectra_dict["data_path"] = self.fake_data_file
        else:
            binned_spectra_dict["data_path"] = self.shapes_files[0]
        binned_spectra_dict["data_format"] = "R counts"
        binned_spectra_dict["data_var_names"] = ["%s_%s"%(self.data_set_names[0],
                                                          self.dim_name)]
        binned_spectra_dict["parameters"] = []
        for i_param in range(self.num_params):
            binned_spectra_dict["parameters"].append(
                {
                    "name": self.param_names[i_param],
                    "shape_path": self.shapes_files[0],
                    "shape_format": "R counts",
                    "shape_var_names": ["%s_%i_%s"%
                                        (self.param_names[i_param], 0,
                                         self.dim_name)],
                    "distribution_path": self.morpho_output_dir+"/"+\
                                         self.morpho_output_file+".root",
                    "distribution_format": "root values",
                    "distribution_tree": self.morpho_output_tree,
                    "distribution_branches": ["rate_"+self.param_names[i_param]]
                }
            )
        plotting.append(binned_spectra_dict)

        for i_data in range(self.num_data_sets):
            binned_spectra_dict = UnsortableOrderedDict()
            binned_spectra_dict["module_name"] = "binned_spectra"
            binned_spectra_dict["method_name"] = "reconstructed_spectrum"
            binned_spectra_dict["output_dir"] = self.plots_output_dir
            binned_spectra_dict["individual_param_output_dir"] = self.plots_output_dir +\
                                                                 "/param_reconstructions"
            binned_spectra_dict["output_path_prefix"] =  self.data_set_names[i_data]+"_"
            binned_spectra_dict["output_format"] = self.plots_output_format
            binned_spectra_dict["store_param_dists"] = False
            binned_spectra_dict["store_param_fractions"] = True
            binned_spectra_dict["store_param_dists_dir"] = self.plots_output_dir+"/param_dists"
            binned_spectra_dict["make_individual_spectra"] = True
            binned_spectra_dict["make_stacked_spectra"] = True
            binned_spectra_dict["make_unstacked_spectra"] = True
            binned_spectra_dict["make_reconstruction_plot"] = True
            binned_spectra_dict["make_diff_plot"] = True
            binned_spectra_dict["make_residual_pull_plot"] = True
            binned_spectra_dict["make_data_model_ratio_plot"] = True
            binned_spectra_dict["make_chi2_vs_dof_plot"] = True
            binned_spectra_dict["make_data_plot"] = True
            binned_spectra_dict["make_param_dist_plots"] = False
            binned_spectra_dict["binning_file"] = self.binning_files[i_data]
            binned_spectra_dict["binning_file_format"] = "R"
            binned_spectra_dict["binning_file_variable"] = self.dim_name
            binned_spectra_dict["divide_by_bin_width"] = True
            binned_spectra_dict["xlabel"] = "Energy (keV)"
            binned_spectra_dict["ylabel"] = "Counts/keV"
            binned_spectra_dict["ylog"] = True
            binned_spectra_dict["title_prefix"] = self.data_set_names[i_data] + " "
            if self.generate_fake_data:
                binned_spectra_dict["data_path"] = self.fake_data_file
            else:
                binned_spectra_dict["data_path"] = self.shapes_files[i_data]
            binned_spectra_dict["data_format"] = "R counts"
            binned_spectra_dict["data_var_names"] = ["%s_%s"%(self.data_set_names[i_data],
                                                              self.dim_name)]
            binned_spectra_dict["parameters"] = []
            for i_param in range(self.num_params):
                binned_spectra_dict["parameters"].append(
                    {
                        "name": self.param_names[i_param],
                        "shape_path": self.shapes_files[i_data],
                        "shape_format": "R counts",
                        "shape_var_names": ["%s_%i_%s"%
                                            (self.param_names[i_param],i_data,
                                             self.dim_name)],
                        "distribution_path": self.morpho_output_dir+"/"+\
                                             self.morpho_output_file+".root",
                        "distribution_format": "root values",
                        "distribution_tree": self.morpho_output_tree,
                        "distribution_branches": ["rate_"+self.param_names[i_param]]
                    }
                )
            plotting.append(binned_spectra_dict)

        '''# Add correlation plots
        apost_dict = UnsortableOrderedDict()
        apost_dict["method_name"] = "aposteriori_distribution"
        apost_dict["module_name"] = "histo"
        apost_dict["input_file_name"] = self.morpho_output_dir+"/"+\
                                        self.morpho_output_file+".root"
        apost_dict["input_tree"] = self.morpho_output_tree
        apost_dict["root_plot_option"] = "cont"
        apost_dict["output_path"] =  self.plots_output_dir
        apost_dict["title"] = "posterior_distros"
        apost_dict["output_format"] = "pdf"
        apost_dict["output_width"] = 12000
        apost_dict["output_height"] = 11000
        apost_dict["data"] = []
        for p_name in self.param_names:
            apost_dict["data"].append("rate_"+p_name)
        plotting.append(apost_dict)'''

        corr_factors_dict = UnsortableOrderedDict()
        corr_factors_dict["module_name"] = "histo"
        corr_factors_dict["method_name"] = "correlation_factors"
        corr_factors_dict["input_file_name"] = self.morpho_output_dir+"/"+\
                                                self.morpho_output_file+".root"
        corr_factors_dict["input_tree"] = self.morpho_output_tree
        corr_factors_dict["output_path"] = self.plots_output_dir
        corr_factors_dict["title"] = "correlation_factors"
        corr_factors_dict["output_format"] = "pdf"
        corr_factors_dict["output_width"] = 12000
        corr_factors_dict["output_height"] = 12000
        corr_factors_dict["data"] = []
        for p_name in self.param_names:
            corr_factors_dict["data"].append("rate_"+p_name)
        plotting.append(corr_factors_dict)
            
        for plot_dict in self.which_plot:
            plotting.append(plot_dict)

        which_plot = UnsortableOrderedDict()
        which_plot["which_plot"] = plotting
        morpho_config_dict["plot"] = which_plot
            

        return yaml.dump(morpho_config_dict, default_flow_style=False)

    def test_morpho_config(self, config_dict):
        """Tries to configure Morpho modules with the given config

        Calls Configure(processor_dict) for each relevant processor,
        and prints any errors that are returned.

        Args:
            config_dict: String containing the yaml dictionary used to
                configure Morpho

        Returns:
            None: Prints all output and any thrown exceptions from
            the tested Morpho processors
        """
        return

    def get_stan_model(self):
        """Build the Stan model

        Args:
            None: Config() must have been run

        Returns:
            str: String containing the stan model code
        """
        d_name = self.dim_name
        model = "/* CUORE analysis Stan model\n"+\
                "* Automatically generated by BinnnedConfigBuilder object\n*/\n\n"+\
                "functions {\n\n}\n\n"

        model += "data {\n\n"+\
                "  int nBins_%s;\n\n"%d_name+\
                "  // Fake Data\n"
        for i_data in range(self.num_data_sets):
            model += "  int %s_%s[nBins_%s];\n"%(self.data_set_names[i_data], d_name, d_name)
        model += "\n"
        for i_data in range(self.num_data_sets):
            model += "  // %s Shapes \n"%self.data_set_names[i_data]
            for p_name in self.param_names:
                model += "  vector[nBins_%s] %s_%i_%s;\n"%(d_name,p_name,i_data,d_name)
            model += "\n"
        model += "}\n\n"

        model += "parameters {\n\n"
        for i_param,p_name in enumerate(self.param_names):
            model += "  real<lower=%s, upper=%s> rate_%s;\n"%\
                     (self.param_lower_bounds[i_param],
                      self.param_upper_bounds[i_param],
                      self.param_names[i_param])
        model += "\n}\n\n"

        model += "transformed parameters {\n\n"
        for i_data in range(self.num_data_sets):
            model += "  real n_counts_recon_%i[nBins_%s];\n"%(i_data,d_name)
        model += "\n"
        for i_data in range(self.num_data_sets):
            model += "  for(i in 1:nBins_%s){\n"%d_name+\
                     "    n_counts_recon_%i[i] =\n"%i_data
            for i_param, p_name in enumerate(self.param_names):
                model += "      rate_%s*%s_%i_%s[i]+\n"%(p_name,p_name,i_data,d_name)
            model = model[:-2] + ";\n"
            model += "  }\n\n"
        model += "}\n\n"

        model += "model {\n\n"
        # Add Priors
        for i_param, p_name in enumerate(self.param_names):
            prior = self.param_priors[i_param]
            if not (prior is None or
                    prior=="None" or
                    prior==""):
                model += "  rate_%s ~ %s;\n"%(p_name, prior)
        model += "\n"
        # Update Likelihood
        model += "  for(i in 1:nBins_%s){\n"%d_name
        for i_data in range(self.num_data_sets):
            model += "    if(n_counts_recon_%i[i]>0){\n"%(i_data)+\
                     "      target += poisson_lpmf(%s_%s[i] | n_counts_recon_%i[i]);\n"%\
                     (self.data_set_names[i_data], d_name,i_data)+\
                     "    }\n"
        model += "  }\n\n}\n"

        return model

    def Run(self):
        """Create the Morpho config and Stan model"""
        if self.generate_morpho_config:
            logger.info("Building Morpho config dictionary")
            morpho_dict = self.get_morpho_config()
            path = os.path.dirname(self.morpho_config_output_path)
            try:
                os.makedirs(path)
            except OSError:
                if not os.path.isdir(path):
                    raise
            morpho_file = open(self.morpho_config_output_path, 'w')
            morpho_file.write(morpho_dict)
            logger.info("Morpho config saved to %s"
                        %self.morpho_config_output_path)
        
            logger.info("Testing Morpho config dictionary")
            self.test_morpho_config(morpho_dict)
            logger.info("Done testing Morpho config")
        else:
            logger.info("Skipping building Morpho config dictionary")

        if self.generate_stan_model:
            logger.info("Building Stan model")
            stan_model = self.get_stan_model()
            path = os.path.dirname(self.stan_model_output_path)
            try:
                os.makedirs(path)
            except OSError:
                if not os.path.isdir(path):
                    raise
            stan_file = open(self.stan_model_output_path, 'w')
            stan_file.write(stan_model)
            stan_file.close()
            logger.info("Stan model saved to %s"
                        %self.stan_model_output_path)
        else:
            logger.info("Skipping building Stan model")

        return

def generate_binned_config(param_dict):
    """Generate Morpho config and Stan model

    Creates a BinnedConfigBuilder object, configures with
    param_dict, and then runs.

    Args:
        param_dict: dictionary used to configure the
            BinnedConfigBuilder object

    Returns:
        None: Generates the config file and Stan model
    """
    proc = BinnedConfigBuilder("config_builder")
    proc.Configure(param_dict)
    proc.Run()
    return
    
if __name__=='__main__':
    p = ArgumentParser(description="Create Morpho config and Stan model files")
    p.add_argument('-c','--config',
                   metavar='<configuration file>',
                   help='Full path to file used to specify configuration',
                   required=True)
    p.add_argument('-v', '--verbosity', default='DEBUG',
                   metavar='<verbosity>',
                   help="Specify verbosity of the logger, with options DEBUG, INFO, WARNING, ERROR, or CRITICAL (Default: DEBUG)",
                   choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'],
                   required=False)
    args = p.parse_args()
    logger.setLevel(args.verbosity)
    with open(args.config, 'r') as cfile:
        param_dict = yaml.load(cfile)
    generate_binned_config(param_dict)



