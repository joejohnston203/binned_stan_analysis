#======================================================
# binned_fake_data.py
#
# Author: J. Johnston
# Date: Mar. 1, 2018
#
# Generate fake data for a binned analysis
#=======================================================

"""Generate fake data for a binned analysis

Classes:
  - BinnedFakeDataProcessor: Class to generate fake data

Functions:
  - binned_fake_data: Generate and store fake data

ToDo:
  - Allow for hierarchical parameters that fluctuate separately in each
    bin, or in each bin along a dimension.
"""

import logging
logger = logging.getLogger(__name__)

import numpy as np

from morpho.utilities.reader import read_param
from morpho.utilities.file_reader import *
from morpho.utilities.file_writer import *

# Implement as a processor in order to prepare for integration with morpho 2
class BinnedFakeDataProcessor:
    """Generate and store fake data

    params is constructed from preprocessing.which_pp in the yaml
    file used to configure morpho.

    Args:
        module_name: Must be "binned_fake_data" to specify this module
        method_name: Must be "binned_fake_data" to specify this 
            method
        generate_fake_data: Boolean specifying whether the method
            should run (Default=True)

        output_dir: String specifying output directory (required)
        output_path_prefix: String specifying output filename
            prefix (Default="fake_data")
        output_format: Format for output data. Options are "text",
            "root", or "R". (Default "R")
        output_variable: Variable used to store the data. For root,
            should be a length 2 list specifying a tree and branch
            name. For R, should be a string specifying the variable
            name. (Default="fake_data")

        poisson_fluctuate: Whether the total number of counts in each
            bin should be poisson fluctuated. (Default=True)

        debug_dir: String specifying directory to store debug text and
            plots. Set to "" to not output debug info. (Default="")

        parameters: List of dictionaries containing the parameters used
            to generate fake data, and the magnitude used to multiply
            each parameter. Each should contain the following keys:
          - name: Name of the parameter (Optional, used for debug output)
          - magnitude: Float used to multiply the shape (required)
          - gaussian_fluctuate: Gaussian fluctuate the magnitude
                (Default=False)
          - gaussian_sigma: Sigma for the fluctuation (Default=0.0)
          - shapes: List of dictionaries used to obtain the shapes.
            numpy.outer will be used to combine the shapes into an n-dim
            array. The dimensionality of all shapes must agree, or else
            fake data cannot be generated. (required). Each dictionary
            must contain the following keys:
              - path: Path to the file containing the shape (required)
              - format: Format of the file. Options are "text", "root",
                and "R" (required)
              - variables: Variables used to access the shape. For
                "text", pass a string specifying elements ("1,:" specifies
                second row, ":,1" specifies second column, ":" or ":,:"
                specifies all).
                For "root", pass a length 2 list,
                ["tree_name","branch_name"].
                For "R", pass a string giving the variable name.

    Returns:
        None: Run() generates and stores fake data to file.
    """
    def __init__(self, name, *args, **kwargs):
        self.__name = name
        return

    def Configure(self, params):
        self.generate_fake_data = read_param(params,
                                             'generate_fake_data',
                                             True)
        if not self.generate_fake_data:
            return

        self.output_dir = read_param(params, 'output_dir', 'required')
        self.output_path_prefix = read_param(params,
                                             'output_path_prefix',
                                             "fake_data")
        self.output_format = read_param(params, 'output_format', 'R')
        if not (self.output_format=="text" or
                self.output_format=="root" or
                self.output_format=="R"):
            logger.warn("Invalid output format %s"%self.output_format)
            logger.warn("Using R")
            self.output_format = "R"
        self.output_variable = read_param(params, 'output_variable',
                                          'fake_data')

        self.poisson_fluctuate = read_param(params, 'poisson_fluctuate',
                                            True)

        self.debug_dir = read_param(params, 'debug_dir', "")
        self.save_debug = (not self.debug_dir=="")

        self.parameters = read_param(params, 'parameters', 'required')

        # Access and store shapes
        self.magnitudes = []
        self.gaussian_fluctuate = []
        self.gaussian_sigma = []
        self.param_shapes = []
        num_shapes = len(read_param(self.parameters[0],
                                    'shapes', 'required'))
        self.dimension_nbins = [None]*num_shapes
        for i,p in enumerate(self.parameters):
            self.param_shapes.append([])
            p["name"] = read_param(p, 'name', "param_%i"%i)
            logger.info("Getting shapes for parameter %s"%p["name"])
            self.magnitudes.append(read_param(p, 'magnitude', 'required'))
            self.gaussian_fluctuate.append(read_param(p,
                                                      'gaussian_fluctuate',
                                                      False))
            self.gaussian_sigma.append(read_param(p,
                                                  'gaussian_sigma',
                                                  0.0))
            shapes = read_param(p, 'shapes', 'required')
            if not len(shapes)==num_shapes:
                logger.error("Parameter %s has %i dimensions, but should have %i"
                             %(p["name"], len(shapes), num_shapes))
                logger.error("Not generating fake data")
                self.generate_fake_data = False
            for j,s in enumerate(shapes):
                path = read_param(s, 'path', 'required')
                file_format = read_param(s, 'format', 'required')
                variables = read_param(s, 'variables', None)
                try:
                    self.param_shapes[i].append(
                        get_variable_from_file(path, file_format, variables))
                    if self.dimension_nbins[j] is None:
                        self.dimension_nbins[j] = len(self.param_shapes[i][j])
                    elif not (len(self.param_shapes[i][j])==
                              self.dimension_nbins[j]):
                        raise ValueError("Dimension mismatch error in shape %s"%j)
                except Exception as e:
                    logger.error("Got error when accessing shape: %s"%e)
                    logger.error("Shape path: %s"%path)
                    logger.error("Shape format: %s"%file_format)
                    logger.error("Variables: %s"%variables)
                    logger.error("Not generating fake data")
                    self.generate_fake_data = False
                    return
        return

    def Run(self):
        if not self.generate_fake_data:
            logger.info("Not generating fake data because generate_fake_data==False")
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        data = np.zeros(tuple(self.dimension_nbins))
        for i in range(len(self.magnitudes)):
            curr_shape = reduce(np.multiply.outer, self.param_shapes[i])
            if self.gaussian_fluctuate[i]:
                magnitude = np.random.normal(self.magnitudes[i],
                                             self.gaussian_sigma[i])
            else:
                magnitude = self.magnitudes[i]
            data += magnitude*curr_shape

        if self.poisson_fluctuate:
            data = np.random.poisson(data)

        file_ending = {"text":".txt", "root":".root", "R":".out"}
        output_path = self.output_dir + "/" + self.output_path_prefix + \
                      file_ending[self.output_format]
        write_variable_to_file(output_path,
                               self.output_format,
                               data,
                               self.output_variable,
                               False)
        logger.info("Fake data stored at %s"%output_path)
        return

def binned_fake_data(param_dict):
    """Generates and stores fake data

    Creates a BinnedFakeDataProcessor object, configures with
    param_dict, and then runs.

    Args:
        param_dict: dicitonary used to configure the
            BinnedFakeDataProcessor object

    Returns:
        None: Generates and stores shapes to file
    """
    proc = BinnedFakeDataProcessor("fake_data")
    proc.Configure(param_dict)
    proc.Run()
    return
