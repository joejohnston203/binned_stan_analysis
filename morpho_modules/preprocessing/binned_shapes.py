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
  - binning_algorithm_default: Determine binning from data
  - generate_shapes: Create and store shapes

ToDo:
  - Add the option to output debugging plots
  - Implement the binning algorithm
"""

import logging
logger = logging.getLogger(__name__)

import numpy as np
import numbers

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

        debug_dir: String specifying directory to store debug text
            and plots. Set to "" to not output debug info. (Default="")

        dimensions: List of dictionaries, with one dictionary for each
            dimension of the analysis (eg "energy" and "time"). (required).
            Each dictionary must contain the following keys:
          - name: Name of the dimension
          - binning_type: How binning should be determined. Options are
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
          - regenerate: Boolean specifying whether this shape should be
            generated.
          - shapes: List of dictionaries, with one dictionary for each
            dimension defined above. Each dictionary should contain
            the following keys:
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

    Returns:
        None: Run() stores the binning and generated shapes to file. All
        outputs are stored in the given output_dir. 
        If output_format=="text", then the binning is stored in text files
        output_path_prefix+"binning_"+name+".txt" and the 
        shapes are stored in
        output_path_prefix+param_name+"_"+dimension_name+".txt",
        where name is the dimension or parameter name.
        If output_format=="root", then the binnings are stored in
        output_path_prefix+"binnings.root", with tree and branch equal
        to the name of the current dimension. Shapes are stored in
        output_path_prefix+"shapes.root", with tree and branch name
        equal to the param_name+"_"+dimension_name.
        If output_format=="R", then the binnings are in the file
        output_path_prefix+"binnings.out" with variable name given by
        the dimension name, and the shapes are in the file
        output_path_prefix+"shapes.out", with the variable names given by
        param_name+"_"+dimension_name.
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

        self.debug_dir = read_param(params, 'debug_dir', "")
        self.save_debug = (not self.debug_dir=="")

        self.dimensions = read_param(params, 'dimensions', 'required')
        self.bin_dicts = []
        for d in self.dimensions:
            self.bin_dicts.append(read_param(d, 'binning', 'required'))
        self.parameters = read_param(params, 'parameters', 'required')

        for p in self.parameters:
            read_param(p, 'name', 'required')
            for s in read_param(p, 'shapes', 'required'):
                read_param(s, 'path', 'required')
                s["renormalize"] = read_param(s, 'renormalize', False)
                s["multiply_shape"] = read_param(s, 'multiply_shape', 1.0)
                s["number_save_type"] = read_param(s, 'number_save_type', 'float64')
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
                    s["variables"]["cut"] = \
                        read_param(s, 'cut', "")
                if "R" in param_format:
                    s["variables"]["variable_names"] = \
                        read_param(s, 'var_names', 'required')
                if "python" in param_format:
                    s["variables"]["path"] = s["path"]
                    s["variables"]["module"] = \
                        read_param(s, 'module', 'required')
                    s["variables"]["method_name"] = \
                        read_param(s, 'function', 'required')
                    s["variables"]["method_options"] = \
                        read_param(s, 'function_options', {})

        return

    def Run(self):
        if not self.generate_shapes:
            logger.info("Not generating shapes because generate_shapes==False")
            return
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        logger.info("Generating and storing binnings")
        self.binnings = list()
        for i_dim, d in enumerate(self.dimensions):
            # Determine binning
            bin_dict = self.bin_dicts[i_dim]
            bin_min = read_param(bin_dict, 'min', None)
            bin_max = read_param(bin_dict, 'max', None)
            bin_min_width = read_param(bin_dict, 'min_bin_width', None)

            bin_include_peaks = (not bin_min is None and
                                 not bin_max is None and
                                 read_param(bin_dict, 'include_peaks', False))
            if bin_include_peaks:
                try:
                    bin_peak_means = \
                        get_variable_from_file(
                            read_param(bin_dict, 'peak_means_path', None),
                            read_param(bin_dict, 'peak_means_format', None),
                            read_param(bin_dict, 'peak_means_variables', None))
                    bin_peak_widths = \
                        get_variable_from_file(
                            read_param(bin_dict, 'peak_widths_path', None),
                            read_param(bin_dict, 'peak_widths_format', None),
                            read_param(bin_dict, 'peak_widths_variables', None))
                except Exception as e:
                    logger.info("Faled to read peaks with exception %s"%e)
                    bin_peak_means = None
                    bin_peak_widths = None
                if (bin_peak_means is None or
                    bin_peak_widths is None):
                    logger.info("Not including peaks because means and/or widths were not provided")
                    bin_include_peaks = False
            if not bin_include_peaks:
                bin_peak_means = []
                bin_peak_widths = []

            bin_merge_low_stats = read_param(bin_dict, 'merge_low_stats', False)
            if bin_merge_low_stats:
                bin_min_counts = read_param(bin_dict, 'min_bin_counts', 'required')
                data_path = read_param(d, 'binning_data_path', 'required')
                data_format = read_param(d, 'binning_data_format', 'required')
                data_variables = read_param(d, 'binning_data_variables', None)
                data = get_variable_from_file(data_path, data_format,
                                              data_variables)
            rebin_regions = []
            for rebin_region_dict in read_param(bin_dict, 'rebin_regions', []):
                try:
                    curr_region = \
                        get_variable_from_file(
                            read_param(rebin_region_dict, 'path', None),
                            read_param(rebin_region_dict, 'format', None),
                            read_param(rebin_region_dict, 'variables', None))
                except Exception as e:
                    logger.info("Failed to read rebin_region with exception %s"%e)
                    curr_region = None
                if not curr_region is None:
                    rebin_regions.append(curr_region)

            binning = []
            if (isinstance(bin_min, numbers.Number) and
                isinstance(bin_max, numbers.Number)):
                binning = [bin_min]
                for i_peak in range(min(len(bin_peak_means),
                                        len(bin_peak_widths))):
                    last_bin = binning[-1]
                    mean = bin_peak_means[i_peak]
                    width = bin_peak_widths[i_peak]
                    peak_min = mean-width
                    peak_max = mean+width

                    if peak_max <= bin_min:
                        continue
                    elif peak_min <= bin_min:
                        binning.append(peak_max)
                        continue

                    if peak_min - last_bin < 0.:
                        binning[-1] = peak_max
                    elif peak_min - last_bin < bin_min_width:
                        if last_bin != bin_min:
                            binning[-1] = 0.5*(last_bin+peak_min)
                        binning.append(peak_max)
                    else:
                        n_divisions = int((peak_min-last_bin)//bin_min_width)
                        step = (peak_min-last_bin)/float(n_divisions)
                        for j in range(1, n_divisions):
                            binning.append(last_bin + j*step)
                        binning.append(peak_min)
                        binning.append(peak_max)
                last_bin = binning[-1]
                if last_bin < bin_max:
                    n_divisions = int((bin_max - last_bin)//bin_min_width)
                    if(n_divisions > 0):
                        step = (bin_max-last_bin)/float(n_divisions)
                        for j in range(1, n_divisions):
                            binning.append(last_bin + j*step)
                    binning.append(bin_max)
            if len(binning)>0:
                binning.sort()

            if bin_merge_low_stats:
                print("Merge low stats here- implement this")

            for rebin_list in rebin_regions:
                rebin_list.sort()
                rebin_min = min(rebin_list)
                rebin_max = max(rebin_list)
                i_bin = 0
                while(i_bin<len(binning) and
                      binning[i_bin]<rebin_min):
                    i_bin += 1
                while(i_bin<len(binning) and
                      binning[i_bin]<rebin_max):
                    binning.pop(i_bin)
                for new_bin in rebin_list:
                    binning.insert(i_bin, new_bin)
                    i_bin += 1
            self.binnings.append(np.array(binning))

            # Store Binning
            d_name = read_param(d, 'name', 'required')
            if self.output_format=="text":
                binning_output_path = self.output_dir + "/" + \
                                      self.output_path_prefix + \
                                      "binning_" + \
                                      d_name + ".txt"
                binning_var_name = None
            elif self.output_format=="root":
                binning_output_path = self.output_dir + "/" + \
                                      self.output_path_prefix + \
                                      "binnings.root"
                binning_var_name = [d_name, d_name]
            else:
                binning_output_path = self.output_dir + "/" + \
                                      self.output_path_prefix + \
                                      "binnings.out"
                binning_var_name = d_name
            write_variable_to_file(binning_output_path, self.output_format,
                                   self.binnings[-1],
                                   binning_var_name, False)
            logger.info("\tBinning for dimension %s stored at %s"%
                        (d_name, binning_output_path))
            binning_n_bins = len(self.binnings[-1])-1
            write_variable_to_file(self.output_dir + "/" + \
                                   self.output_path_prefix + "binnings.out",
                                   "R", binning_n_bins, "nBins_%s"%d_name)


        logger.info("Generating and storing parameter shapes")
        self.param_shapes = []
        for p in self.parameters:
            p_name = read_param(p, 'name', 'required')
            if not read_param(p, 'regenerate', True):
                logger.info("Skipping shape %s"%p_name)
                continue

            shapes = []
            for i,s in enumerate(read_param(p, 'shapes', 'required')):
                curr_shape = get_histo_shape_from_file(self.binnings[i],
                                                        read_param(s, 'path', 'required'),
                                                        read_param(s, 'format', 'required'),
                                                        read_param(s, 'variables', 'required'))[0]
                if s["renormalize"]:
                    curr_shape = curr_shape.astype('float64')/float(np.sum(curr_shape))
                curr_shape = curr_shape.astype('float64')*float(s["multiply_shape"])
                shapes.append(curr_shape.astype(s["number_save_type"]))
            self.param_shapes.append(shapes)

            if self.output_format=="text":
                for i in range(len(shapes)):
                    shape_output_path = self.output_dir + "/" + \
                                        self.output_path_prefix + \
                                        p_name + "_" + \
                                        self.dimensions[i]["name"] + ".txt"
                    shape_var_name = None
                    write_variable_to_file(shape_output_path, self.output_format,
                                           shapes[i],
                                           shape_var_name, False)
                    logger.info("\tShape for parameter %s, dimension %s stored at %s"%
                                (p_name, self.dimensions[i]["name"], shape_output_path))
            elif self.output_format=="root":
                shape_output_path = self.output_dir + "/" + \
                                    self.output_path_prefix + \
                                    "shapes.root"
                for i in range(len(shapes)):
                    tree_name = p_name + "_" + self.dimensions[i]["name"]
                    branch_names = [tree_name]
                    arrs = [shapes[i]]
                    write_root_branches(shape_output_path, tree_name,
                                        branch_names, arrs, False)
                logger.info("\tShapes for parameter %s stored at %s"%
                            (p_name, shape_output_path))
            else:
                shape_output_path = self.output_dir + "/" + \
                                    self.output_path_prefix + \
                                    "shapes.out"
                for i in range(len(shapes)):
                    shape_var_name = p_name + "_" + \
                                     self.dimensions[i]["name"]
                    write_variable_to_file(shape_output_path, self.output_format,
                                           shapes[i], shape_var_name, False)
                logger.info("\tShapes for parameter %s stored at %s"%
                            (p_name, shape_output_path))
        return

def binning_algorithm_default(lb, ub, data, args):
    """Algorithm to determine binning given data

    Args:
        lb, ub: Lower bound and upper bound for the binning
        data: List of data point to histogram when determining
            binning
        args: Other arguments (will be expanded later)
    """
    logger.warn("binning algorithm not yet implemented")
    logger.warn("Using 20 uniform bins between %s, %s"%(lb,ub))
    return np.linspace(lb, ub, 21, endpoint=True)

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
