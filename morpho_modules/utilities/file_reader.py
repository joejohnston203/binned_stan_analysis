#======================================================
# binned_spectra.py
#
# Author: J. Johnston
# Date: Mar. 1, 2018
#
# Read data from different file types
#=======================================================

import logging
logger = logging.getLogger(__name__)

import numpy as np
import pystan

try:
    import ROOT
except ImportError as e:
    logger.info("ROOT could not be imported, error: %s"%e)
    logger.info("Continuing without ROOT")

def read_txt_column(file_path, column):
    return []

def read_root_branch(file_path, tree, branch):
    return []

def read_R_variable(file_path, var_name):
    return []

def evaluate_python_fcn(xvals, module_name, fcn_name):
    return []

def get_variable_from_file(path, file_format, variable=None):
    """Get a variable from file

    Args:
        path: Path to the file
        file_format: Format of the file. Options are "text",
            "root", "R", or "python"
        variable: Variable to access. The form depends on the file
            format:
          - "text": The file will be loaded with np.loadtxt. Set
            variable=None to load the entire array, or pass an int
            specifying the column to load.
          - "root": variable is a length 2 list of str specifying the
            tree name and branch name.
          - "R": variable is a string specifying the variable name
          - "python": variable is a length 2 list specifying the
            the path to the module and the function name

    Returns:
        Depends on the file format:
          - "text", "root": np.array containing the specified array
          - "R": Returns the specified variable (can be any type)
          - "python": python function
    """
    if file_format=="text":
        print("reading text not yet implemented")
    elif file_format=="root":
        print("Reading root not yet implemented")
    elif file_format=="R":
        r_dict = pystan.misc.read_rdump(path)
        res_variable = r_dict[variable]
    elif file_format=="python":
        print("reading python not yet implemented")
    else:
        logger.warn("Invalid file_format. Returning None")
        res_variable = None
    return res_variable
    

def get_histo_shape_from_file(binning, path, file_format, variables):
    """Get a histogram shape from file

    Args:
        binning: List giving the bin edges for the histogram
        path: String giving path to the data.
        file_format: Format used to save the data. Currently supported
            options are "text values", "text counts", "text function",
            "root values", "root counts", "root function", "R values",
            "R counts", "R function", "python function".
            "text", "root", "R", or "python" specifies the filetype.
            "values" specifies that there will be a 1D list of values to
            be histogrammed. "counts" specifies that the number of counts
            in each bin will be specified (assuming bins from binning_file
            or n_bins). "function" specifies that there will be two columns
            of points defining the spectrum, which will then be integrated
            in each bin.
        variables: Dictionary containing the variables used to access
            the data. The formatting depends on the file format:
          - "text": variables should containg a key "columns",
            containing a list of columns containing the data. There
            should be one column for "values" or "counts",and two
            columns for "function".
          - "root": variables should contain a key "tree" with
            the name of the tree, and a key "branches" with a list of
            branches (one for "values" or "counts", two for "function")
          - "R": variables should contain a key "variable_names" with
            a list of variable names (one for "values" or "counts", two
            for "function".
          - "python": variables should contain a key "module" and a key
            "method_name" specifying the path to the module to load, and
            the name of the method defining the function.
    Returns:
        A np.array specifying the number of counts in each bin described by
        binning. If the format was a function, then the value for each bin
        is the number of counts per unit x.
    """
    #print("binning: %s"%binning)
    #print("path: %s"%path)
    #print("file_format: %s"%file_format)
    #print("variables: %s"%variables)

    if "text" in file_format:
        print("text not yet implemented")
    elif "root" in file_format:
        print("root not yet implemented")
    elif "R" in file_format:
        if "counts" in file_format or "values" in file_format:
            res_arr = get_variable_from_file(path, "R", variables["variable_names"][0])
    elif "python" in file_format:
        print("python not yet implemented")
    else:
        logger.warn("Invalid file_format given (%s), returning None")
        return None

    # For values, histogram the given array
    if "values" in file_format:
        pass

    # For a spectrum, integrate over each bin
    if "function" in file_format:
        pass

    return np.array(res_arr)
