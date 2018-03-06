#======================================================
# file_writer.py
#
# Author: J. Johnston
# Date: Mar. 1, 2018
#
# Write data to different file types
#=======================================================

import logging
logger = logging.getLogger(__name__)

import os

import numpy as np
import pystan

try:
    import ROOT
except ImportError as e:
    logger.info("ROOT could not be imported, error: %s"%e)
    logger.info("Continuing without ROOT")

from array import array

def write_txt_array(file_path, arr):
    """Write an array to text with np.savetxt

    Only works for a 1D or 2D array. Simply calls np.savetxt.
    This method exists for the sake of symmetry with
    morpho.utilities.file_reader.read_txt_array.

    Args:
        file_path: Path to the output text file
        arr: 1D or 2D numpy array to save

    Returns:
        None: Saves the array to file

    Throws:
        IOError: If the given file cannot be opened
    """
    np.savetxt(file_path, arr)
    return

def write_root_branches(file_path, tree_name, branch_names, arrs,
                        recreate=False):
    """Write a list to a branch in a root file

    The given tree is overwritten if it already exists

    Args:
        file_path: Path to the root file
        tree_name: Name of the tree to store in
        branch_names: List of branches to write
        arrs: List of arrays to save to the branches. All arrays must
            have the same length. Also, len(branch_names)==len(arrs).
        recreate: Boolean specifying whether the root file should
            be overwritten

    Returns:
        None: Saves the arrays to file

    Throws:
        IOError: If the given file cannot be opened
    """
    if not len(branch_names)==len(arrs):
        logger.error("branch names (%s) has different length than arrs (%s)"%
                     (branch_names, arrs))
        logger.error("Not saving to root file")
        return

    if recreate:
        myfile = ROOT.TFile(file_path, "RECREATE")
    else:
        myfile = ROOT.TFile(file_path, "UPDATE")
    tree = ROOT.TTree(tree_name, tree_name)
    tmp_vals = []
    for i in range(len(arrs)):
        tmp_vals.append(array('f',[ 0 ]))
        tree.Branch(branch_names[i], tmp_vals[i], branch_names[i]+'/F')
    for i in range(len(arrs[0])):
        for j in range(len(arrs)):
            if i<len(arrs[j]):
                tmp_vals[j][0] = arrs[j][i]
        tree.Fill()
    tree.Write()
    myfile.Close()
    return

def write_R_variable(file_path, var_name, variable,
                    recreate=False):
    """Write a variable to an R file

    Args:
        file_path: Path to the R file
        var_name: Name to store the variable under
        variable: Variable to store
        recreate: Whether the R file should be overwritten

    Returns:
        None: Stores the variable in an R file

    Throws:
        IOError: If the given file cannot be opened
    """
    if recreate or not os.path.isfile(file_path):
        r_dict = {}
    else:
        try:
            r_dict = pystan.misc.read_rdump(file_path)
        except ValueError as e:
            logger.critical("Could not open R file due to error: %s"%e)
            logger.critical("Output path: %s"%file_path)
            logger.critical("Output variable name: %s"%var_name)
            logger.critical("Output value: %s"%variable)
            logger.critical("Exiting without writing")
            return
    r_dict[var_name] = variable
    pystan.misc.stan_rdump(r_dict, file_path)
    return

def write_variable_to_file(path, file_format, variable,
                           variable_name=None,
                           recreate=False):
    """Write a variable to file

    Args:
        path: Path to the file
        file_format: Format of the file. Options are "text",
            "root", or "R"
        variable: Variable to store. Must be an array for "text" or
            "root".
        variable_name: Variable name used to store the variable. The
            form depends on the file format:
          - "text": No variable_name required
          - "root": variable_name is a length 2 list of str specifying 
            the tree name and branch name.
          - "R": variable_name is a string specifying the variable name
        recreate: Boolean specifying whether the file should be
            recreated. (Only applies to "root" and "R").

    Returns:
        None: Stores the given variable to file

    Throws:
        IOError: If the given file cannot be opened
    """
    if file_format=="text":
        write_txt_array(path, variable)
    elif file_format=="root":
        write_root_branches(path,
                            variable_name[0], [variable_name[1]],
                            [variable], recreate)
    elif file_format=="R":
        write_R_variable(path, variable_name, variable, recreate)
    else:
        logger.warn("Invalid file_format. Not storing variable %s in path %s"%
                    (variable, path))
    return
