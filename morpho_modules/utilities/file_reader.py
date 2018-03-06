#======================================================
# file_reader.py
#
# Author: J. Johnston
# Date: Mar. 1, 2018
#
# Read data from different file types
#=======================================================

"""Read data from files

Functions:
  - read_txt_array: Read an array stored with np.savetxt
  - read_root_branch: Read a branch from a root file
  - histogram_root_branch: Generate histogram from a root branch
  - read_R_variable: Read a variable from a R file
  - evalute_python_fcn: Evaluate a python function for given values
  - get_variable_from_file: Get a variable from various file types
  - get_histo_shape_from_file: Get a histogram from a file
"""

import logging
logger = logging.getLogger(__name__)

import sys
from importlib import import_module

import numpy as np
import pystan

try:
    import ROOT
except ImportError as e:
    logger.info("ROOT could not be imported, error: %s"%e)
    logger.info("Continuing without ROOT")

def read_txt_array(file_path, idx=None):
    """Read an array from a text file created with np.savetxt

    Only works for 1D or 2D array (Such as those saved by np.savetxt)

    Args:
        file_path: Path to the text file
        idx: String specifying elements to return. None returns the
            entire array. "1,2" will return the element at index (1,2).
            "1" or "1,:" will return the second row, and ":,1" will
            return the second column.

    Returns:
        np.array: Array with the selected column

    Throws:
        IOError: If the given file does not exist
    """
    result = np.loadtxt(file_path)

    if idx is None:
        return result
    else:
        indices = idx.split(',')
        if len(indices)==1:
            if indices[0]==":":
                return result
            else:
                return result[int(indices[0])]
        elif len(indices)==2:
            if indices[0]==":":
                if indices[1]==":":
                    return result
                else:
                    return result[:,int(indices[1])]
            else:
                if indices[1]==":":
                    return result[int(indices[0]),:]
                else:
                    return result[int(indices[0]),
                                  int(indices[1])]
        else:
            logger.error("Invalid index %s, returning full array"%idx)
            return result

def read_root_branch(file_path, tree_name, branch_name):
    """Get a branch from a root file

    Args:
        file_path: Path to the root file
        tree_name: Name of the tree to access
        branch_name: Name of branch to return

    Returns:
        list: Containing all elements of the branch

    Throws:
        IOError: If the given file does not exist
    """
    myfile = ROOT.TFile(file_path,"READ")
    tree = myfile.Get(tree_name)
    result = []
    for elt in tree:
        result.append(getattr(elt,branch_name))
    myfile.Close()
    return result

def histogram_root_branch(binning, file_path, tree_name, branch_name):
    """Get and histogram a branch from a root file

    Args:
        binning: List defining bin endges
        file_path: Path to the root file
        tree_name: Name of the tree to access
        branch_name: Name of branch to return

    Returns:
        3-tuple: (np array containing the counts for each bin,
            mean, sigma)

    Throws:
        IOError: If the given file does not exist
    """
    myfile = ROOT.TFile(file_path,"READ")
    tree = myfile.Get(tree_name)
    rh = ROOT.TH1F("rh", "", len(binning)-1, binning)
    tree.Draw(branch_name+">>rh", "", "goff")
    mean = rh.GetMean()
    sigma = rh.GetStdDev()
    bin_contents = np.empty(len(binning)-1)
    for i in range(len(binning)-1):
        bin_contents[i] = rh.GetBinContent(i+1)
    myfile.Close()
    return (bin_contents, mean, sigma)

def read_R_variable(file_path, var_name):
    """Read an array from an R file

    Args:
        file_path: Path to the R file
        var_name: Name of the variable to return

    Returns:
        Variable from R file

    Throws:
        IOError: If the given file does not exist
        KeyError: If the given variable is not in the R file
    """
    r_dict = pystan.misc.read_rdump(file_path)
    return r_dict[var_name]

def evaluate_python_fcn(xvals, path, module_name, fcn_name, options={}):
    sys.path.append(path)
    temp_module = import_module(module_name)
    fcn=getattr(temp_module,fcn_name)
    fcn = np.vectorize(fcn)
    return fcn(np.array(xvals), **options)

def get_variable_from_file(path, file_format, variable=None):
    """Get a variable from file

    Args:
        path: Path to the file
        file_format: Format of the file. Options are "text",
            "root", or "R"
        variable: Variable to access. The form depends on the file
            format:
          - "text": The file will be loaded with np.loadtxt. Pass a
            string specifying elements to return. None returns the
            entire array. "1,2" will return the element at index (1,2).
            "1" or "1,:" will return the second row, and ":,1" will
            return the second column.
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

    Throws:
        IOError: If the given file does not exist
    """
    if file_format=="text":
        res_variable = read_txt_array(path, variable)
    elif file_format=="root":
        res_variable = read_root_branch(path, variable[0], variable[1])
    elif file_format=="R":
        res_variable = read_R_variable(path, variable)
    else:
        logger.warn("Invalid file_format. Returning None")
        res_variable = None
    return res_variable
    

def get_histo_shape_from_file(binning, path, file_format, variables,
                              get_mean_sigma=True):
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
            in each bin. The points defining the spectrum must be denser
            than the binning. The average value in each bin will be
            returned.
        variables: Dictionary containing the variables used to access
            the data. The formatting depends on the file format:
          - "text":  variables should contain a key "columns" specifying
            elements to return. None returns the
            entire array. "1,2" will return the element at index (1,2).
            "1" or "1,:" will return the second row, and ":,1" will
            return the second column.
          - "root": variables should contain a key "tree" with
            the name of the tree, and a key "branches" with a list of
            branches (one for "values" or "counts", two for "function")
          - "R": variables should contain a key "variable_names" with
            a list of variable names (one for "values" or "counts", two
            for "function".
          - "python": variables should contain a key "path", a key "module"
            and a key "method_name" specifying the path to the module to load,
            and the name of the method defining the function. It may also
            contain a key "method_options", which contains a dictionary
            with named arguments to pass the method.
        get_mean_sigma: Whether the mean and sigma should be calculated.
            0 is returned for both if set to False. (Default True).
    Returns:
        A 3-tupe containing (np.array, float64, float64).
        np.array specifies the number of counts in each bin described by
        binning. If the format was a function, then the value for each bin
        is the number of counts per unit x. The first float64 gives the
        mean of the histogram, and the second gives the standard deviation.

    Throws:
        IOError: If the given file does not exist
    """
    if "text" in file_format:
        if "counts" in file_format or "values" in file_format:
            res_arr = read_txt_array(path, variables["columns"])
    elif "root" in file_format:
        if "counts" in file_format or "values" in file_format:
            res_arr = get_variable_from_file(path, "root",
                                             [variables["tree"], variables["branches"][0]])
        elif "function" in file_format:
            x_arr = get_variable_from_file(path, "root",
                                           [variables["tree"], variables["branches"][0]])
            y_arr = get_variable_from_file(path, "root",
                                           [variables["tree"], variables["branches"][1]])
    elif "R" in file_format:
        if "counts" in file_format or "values" in file_format:
            res_arr = get_variable_from_file(path, "R", variables["variable_names"][0])
        elif "function" in file_format:
            x_arr = get_variable_from_file(path, "R", variables["variable_names"][0])
            y_arr = get_variable_from_file(path, "R", variables["variable_names"][1])
    elif "python" in file_format:
        if not "function" in file_format:
            logger.error("file_format %s is invalid, assuming 'python function'")
            file_format = "python function"
        x_arr = np.linspace(binning[0], binning[-1],
                            50*(len(binning)-1), endpoint=False)
        if not "method_options" in variables:
            variables["method_options"] = {}
        y_arr = evaluate_python_fcn(x_arr, variables["path"],
                                    variables["module"], variables["method_name"],
                                    variables["method_options"])
    else:
        logger.error("Invalid file_format given (%s), returning None"%file_format)
        return None

    average = 0.
    sigma = 0.

    # For values, histogram the given array
    if "values" in file_format:
        if(get_mean_sigma):
            average = np.mean(res_arr)
            sigma = np.std(res_arr)
        res_arr = np.histogram(res_arr, binning)[0]

    # For a spectrum, integrate over each bin
    if "function" in file_format:
        if not len(x_arr)==len(y_arr):
            logger.error("x and y arrays have different lengths. Returning None")
            return None
        bin_counts = np.histogram(x_arr, binning)[0]
        bin_total = np.histogram(x_arr, binning, weights=y_arr)[0]
        bin_avg = bin_total.astype('float64')/bin_counts.astype('float64')
        bin_widths = []
        for i in range(len(binning)-1):
            bin_widths.append(binning[i+1]-binning[i])
        res_arr = bin_avg*np.array(bin_widths)

    if "counts" in file_format or "function" in file_format:
        if not len(binning)-1==len(res_arr):
            logger.error("Counts length different from binning. Returning.")
        elif(get_mean_sigma):
            total = 0.
            nelts = 0
            sigma = 0.0
            bin_centers = []
            for i in range(len(binning)-1):
                bin_centers.append(binning[i] +
                                   (binning[i+1]-binning[i])/2.0)
            for i in range(len(res_arr)):
                total += res_arr[i]*bin_centers[i]
                nelts += res_arr[i]
            average = total/float(nelts)
            if nelts>1:
                for i in range(len(res_arr)):
                    sigma += (bin_centers[i]-average)**2*res_arr[i]
                sigma = np.sqrt(sigma/(nelts-1))
            else:
                sigma = float('inf')

    return (np.array(res_arr), average, sigma)
