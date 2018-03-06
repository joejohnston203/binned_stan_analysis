#======================================================
# binned_shapes.py
#
# Author: J. Johnston
# Date: Mar. 1, 2018
#
# Create and store shapes used for a binned analysis
#=======================================================

"""Create and store shapes used for a binned analysis

The Stan model needs the number of expected counts in each bin in
order to calculate the likelihood for each bin. These methods
obtain the expected number of counts from a variety of file
types, then store the shapes for access by the stan model

Classes:
  - GenerateShapesProcessor: Class to create and store shapes

Functions:
  - generate_shapes: Create and store shapes
"""

import logging
logger = logging.getLogger(__name__)

import numpy as np

from morpho.utilities.reader import read_param
from morpho.utilities.file_reader import *
from morpho.utilities.file_writer import *

# Implement as a processor in order to prepare for integration with morpho 2
class GenerateShapesProcessor:
    """Generate and store shapes from info stored in files

    params is constructed from preprocessing.which_pp in the yaml
    file used to configure morpho.

    Args:
        module_name: Must be "binned_shapes" to specify this module
        method_name: Must be "generate_shapes" to specify this 
            method
        generate_shapes: Boolean specifying whether the method should
            run (Default=True)
        
        output_dir: String specifying output directory (required)
        output_path_prefix: String specifying output path prefix
            (required)
        output_format: Format for output shapes. Options are
            "text", "root", or "R".

        dimensions: List of dictionaries, with one dictionary for each
            dimension of the analysis (eg "energy" and "time"). (required).
            Each dictionary must contain the following keys:
          - name: Name of the dimension
          - binning_type: How binning should be determine. Options are
            "file", "uniform", or "algorithm". If file, then a file must
            give the N+1 bin edges defining the bins. If "uniform", then
            N equally spaced bins between a lower bound and upper bound
            will be used. If "algorithm", then a given algorithm will be
            used to determine the binning based on a given dataset.
          - binning_path: Path to file containing binning (required if
            binning_type=="file")
          - binning_format: Format of binning file. Options are "text",
            "root", and "R". (required for "file")
          - binning_variables: Variables used to access the binning. For
            "text", pass a string specifying elements ("1,:" specifies
            second row, ":,1" specifies second column, ":" or ":,:" specifies all).
            For "root", pass a length 2 list, ["tree_name","branch_name"].
            For "R", pass a string giving the variable name.
          - binning_columns, binning_tree, binning_branches, 
            binning_var_names: Info to access the binning data from the
            given file. See shapes settings below. Required for "text".
          - n_bins: Number of bins (required for "uniform")
          - lower_bound, upper_bound: Lower and upper bound for the
            dimension (required for "uniform" or "algorithm")
          - binning_algorithm: Algorithm used to determine the binning
            based on a dataset. Currently the only option is "default".
          - binning_data_path: Path to data used by the algorithm to
            determine binning. The data should be a list of values that will
            be histogrammed. (required for "algorithm").
          - binning_data_format: Format of the data file. Options are "text",
            "root", and "R". (required for "algorithm")
          - binning_data_variables: Variables used to access the data. For
            "text", pass a string specifying elements ("1,:" specifies
            second row, ":,1" specifies second column, ":" or ":,:" specifies all).
            For "root", pass a length 2 list, ["tree_name","branch_name"].
            For "R", pass a string giving the variable name. (required for "algorithm")

        parameters: List dictionaries containing information about the 
            parameters that shapes should be generated for. (required).
            Each dictionary should contain the following keys:
          - name: Name of the parameter (used in saving shapes)
          - regnerate: Boolean specifying whether this shape should be
            generated.
          - shapes: List of dictionaries, with one dictionary for each
            dimension defined above. Each dictionary should contain
            the following keys:
              - path: Path to the file containing the shape
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
                Required if format includes "python".
              - function: Function name. Required if format includes "python".

    Returns:
        None: Run() stores the binning and generated shapes to file. All
        outputs are stored in the given output_dir. 
        If output_format=="text", then the binning is stored in text files
        output_path_prefix+"binning_"+name+".txt" and the 
        shapes are stored in output_path_prefix+name+".txt",
        where name is the dimension or parameter name.
        If output_format=="root", then the binnings are stored in
        output_path_prefix+"binnings.root" and shapes are stored in
        output_path_prefix+"shapes.root". The tree and branch name
        are equal to the name of the current dimension or parameter.
        If output_format=="R", then the binnings are in the file
        output_path_prefix+"binnings.out", and the shapes are in the file
        output_path_prefix+"shapes.out", with the variable names given by
        the dimension and parameter names.
    """
    def __init__(self, name, *args, **kwargs):
        self.__name = name
        return

    def Configure(self, params):
        self.generate_shapes = read_param(params, 'generate_shapes', True)
        if not self.generate_shapes:
            return
        
        self.output_dir = read_param(params, 'output_dir', 'required')
        self.output_path_prefix = read_param(params, 'output_path_prefix', 'required')
        self.output_format = read_param(params, 'output_format', 'R')
        if not (self.output_format=="text" or
                self.output_format=="root" or
                self.output_format=="R"):
            logger.warn("Invalid output format %s"%self.output_format)
            logger.warn("Using R")
            self.output_format = "R"

        self.dimensions = read_param(params, 'dimensions', 'required')
        self.parameters = read_param(params, 'parameters', 'required')

        for d in self.dimensions:
            read_param(d, 'name', 'required')
            binning_type = read_param(d, 'binning_type', 'required')
            if binning_type=="file":
                read_param(d, 'binning_path', 'required')
                read_param(d, 'binning_format', 'required')
                d['binning_variables'] = read_param(d, 'binning_variables', None)
            elif binning_type=="uniform":
                read_param(d, 'n_bins', 'required')
                read_param(d, 'lower_bound', 'required')
                read_param(d, 'upper_bound', 'required')
            elif binning_type=="algorithm":
                read_param(d, 'binning_data_path', 'required')
                read_param(d, 'binning_data_format', 'required')
                d["binning_data_variables"] = read_param(d, 'binning_data_variables', None)
            else:
                logger.error("Invalid binning_type: %s"%binning_type)
                logger.error("Not generating shapes")
                self.generate_shapes = False
                return
            print(d)

        for p in self.parameters:
            read_param(p, 'name', 'required')
            for s in read_param(p, 'shapes', 'required'):
                read_param(s, 'path', 'required')
                param_format = read_param(s, 'format', 'required')
                s["variables"] = {}
                if "text" in param_format:
                    s["variables"]["columns"] = \
                        read_param(s, 'columns', 'required')
                if "root" in param_format:
                    s["variables"]["tree"] = \
                        read_param(s, 'tree', 'required')
                    s["variables"]["branches"] = \
                        read_param(s, 'branches', 'required')
                if "R" in param_format:
                    s["variables"]["variable_names"] = \
                        read_param(s, 'var_names', 'required')
                if "python" in param_format:
                    s["variables"]["module"] = \
                        read_param(s, 'module', 'required')
                    s["variables"]["method_name"] = \
                        read_param(s, 'function', 'required')

        return

    def Run(self):
        if not self.generate_shapes:
            return
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # Set up binning
        self.binning = list()
        for d in self.dimensions:
            binning_type = read_param(d, 'binning_type', 'required')
            if binning_type=="file":
                binning_path = read_param(d, 'binning_path', 'required')
                binning_format = read_param(d, 'binning_format', 'required')
                binning_variables = d["binning_variables"]
                binning = get_variable_from_file(binning_path, binning_format,
                                                 binning_variables)
                self.binning.append(np.array(binning))
            elif binning_type=="uniform":
                n_bins = read_param(d, 'n_bins', 'required')
                lb = read_param(d, 'lower_bound', 'required')
                ub = read_param(d, 'upper_bound', 'required')
                binning = np.linspace(lb, ub, n_bins+1, endpoint=True)
                self.binning.append(binning)
            elif binning_type=="algorithm":
                binning_algorithm = read_param(d, 'binning_algorithm', "default")
                if not binning_algorithm=="default":
                    logger.info("Currently only the default binning algorithm is implemented. Using default.")
                lb = read_param(d, 'lower_bound', 'required')
                ub = read_param(d, 'upper_bound', 'required')
                data_path = read_param(d, 'binning_data_path', 'required')
                data_format = read_param(d, 'binning_data_format', 'required')
                data_variables = read_param(d, 'binning_data_variables', None)
                data = get_variable_from_file(data_path, data_format,
                                              data_variables)
                self.binning.append(binning_algorithm_default(lb, ub, data, None))
            if self.output_format=="text":
                binning_output_path = self.output_dir + "/" + \
                                      self.output_path_prefix + \
                                      "binning_" + \
                                      read_param(d, 'name', 'required') + ".txt"
                binning_var_name = None
            elif self.output_format=="root":
                binning_output_path = self.output_dir + "/" + \
                                      self.output_path_prefix + \
                                      "binnings.root"
                binning_var_name = [read_param(d, 'name', 'required'),
                                    read_param(d, 'name', 'required')]
            else:
                binning_output_path = self.output_dir + "/" + \
                                      self.output_path_prefix + \
                                      "binnings.out"
                binning_var_name = read_param(d, 'name', 'required')
            write_variable_to_file(binning_output_path, self.output_format,
                                   self.binning[-1],
                                   binning_var_name, False)

        for s in self.parameters:
            pass
        return

def binning_algorithm_default(lb, ub, data, args):
    """Algorithm to determine binning given data

    Args:
        lb, ub: Lower bound and upper bound for the binning
        data: List of data point to histogram when determining
            binning
        args: Other arguments (will be expanded later)
    """
    print("binning data = %s"%data)
    print("binning algorithm not yet implemented")
    return []

def generate_shapes(param_dict):
    """Generates and stores binned spectra shapes

    Creates a GenerateShapesProcessor object, configures with
    param_dict, and then runs.

    Args:
        param_dict: dicitonary used to configure the
            GenerateShapesProcessor object

    Returns:
        None: Generates and stores shapes to file
    """
    proc = GenerateShapesProcessor("shapes_generator")
    proc.Configure(param_dict)
    proc.Run()
    return
