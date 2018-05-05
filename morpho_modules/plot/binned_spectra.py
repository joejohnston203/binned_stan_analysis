#======================================================
# binned_spectra.py
#
# Author: J. Johnston
# Date: Mar. 1, 2018
#
# Make plot relevant for an analysis of a binned
# spectrum
#=======================================================

"""Make plots relevant to a binned spectrum analysis

Classes:
  - ReconstructSpectrumProcessor: Class to plot reconstruction

Functions:
  - reconstructed_spectrum: Plot reconstruction with stan parameters

ToDo:
  - Display error bars on the plots
  - Enable specifying a column or row from an R array
  - All plots currently assume 1D. Maybe I can naturally extend to
    multiple dimensions?
  - Output the mean and stddev for each parameter. Possibly calculate
    the rate per region for each parameter. That is, I think if I
    get the sum over the shape for each shape, then divided by it,
    I would get the total expected counts. If I
    then divide by the range of each dimension, I would get the rate
    in counts/keV/day (for example)
"""

import copy
import os

import logging
logger = logging.getLogger(__name__)

import numpy as np

from morpho.utilities.reader import read_param
from morpho.utilities.list_plotter import plot_curves
from morpho.utilities.file_reader import *

# Implement as a processor in order to prepare for integration with morpho 2
class ReconstructSpectrumProcessor:
    """Create a spectrum reconstructed from stan parameters

    params is constructed from plot.which_plot in the yaml
    file used to configure morpho. 

    Args:
        module_name: Must be "binned_spectra" to specify this module
        method_name: Must be "reconstructed_spectrum" to specify this 
            method

        output_dir: String specifying output directory
        output_path_prefix: String specifying output path prefix

        plot_data: Boolean specifying whether data should be
            included on plots (Default=True)
        make_individual_spectra: Boolean specifying whether plots
            should be saved with each spectrum plotted separately
            (Default=True)
        make_stacked_spectra: Boolean specifying whether a plot should
            be saved with all spectra on the same plot (Default=True)
        make_unstacked_spectra: Boolean specifying whether a plot should
            be saved with all spectra on the same plot, stacked such that
            the summed spectra should approximately add to the data.
            (Default=True)
        make_data_plot: Whether a plot should be made containing only
            the data. (Default=True)

        binning_file: Path to text file containing edges of all bins.
            N+1 edges must be specified 
        n_bins: Number of bins to use if the binning file is not
            provided, or if the file does not exist. Regular binnning
            between xmin and xmax is assumed. (Default=50)
        divide_by_bin_width: Boolean specifying whether the counts in
            each bin should be divided by the bin width (Default=True)

        xmin: Minimum x value to plot (Required only if n_bins is used)
        xmax: Maximum x value to plot (Required only if n_bins is used)
        xlabel: String with x axis label (Default="x")
        xlog: Boolean specifying if x is log scaled (Default=False)
        ybounds: Two element list with y bounds (Default="auto")
        ylabel: String with y axis label (Default="Counts Per Unit x")
        ylog: Boolean specifying if y is log scaled (Default=False)
        title_prefix: String with the prefix for plot titles

        data_path: String giving path to the data.
        data_format: Format used to save the data. Currently supported
            options are "text values", "text counts", "text function",
            "root values", "root counts", "root function", "R values",
            "R counts", "R function", "python function".
            "text", "root", "R", or "python" specifies the filetype.
            "values" specifies that there will be a 1D list of values to
            be histogrammed. "counts" specifies that the number of counts
            in each bin will be specified (assuming bins from binning_file
            or n_bins). "function" specifies that there will be two columns
            of points defining the spectrum.
        data_columns: List of column(s) including the relevant data. One column
            number should be included for "values" or "counts", two for "function".
            Required  if format includes "text".
        data_tree: Root tree name used to access data. Required if format
            includes "root".
        data_branches: list of root branch name(s) used to access data. 
            (See data_columns). Required if format includes "root".
        data_var_names: Variable names in the "R" file. (See data_columns).
            Required if format includes "R".
        data_module: Module where the python function will be found.
            Required if format includes "python". In this case, data_path
            specifies the the path to the directory where the module is
            located.
        data_function: Function name. Required if format includes "python".
        data_function_options: Dictionary with any named arguments to pass
            to the function. (Default={})

        parameters: Array of dictionaries containing the following keys
          - shape_path: Path to the file where the shape is stored
          - shape_format: Format used to save the shape (see data_format)
          - shape_columns: Columns with shape (see data_columns)
          - shape_tree: Root tree used to access shape (see data_tree)
          - shape_branches: Root branches used to access shape 
            (see data_branches)
          - shape_var_names: R variable names (see data_var_names)
          - shape_module: Module name (see data_module)
          - shape_function: Function name (see data_function)
          - shape_function_options: Options (see data_function_options)
          - distribution_path: Path to file with the distribution
            of the parameter fromm the sampler
          - distribution_format: Format of file with the distribution.
            Note that "function" is not implemented for distributions.
          - distribution_columns: Columns with the distribution
          - distribution_tree: Root tree used to access distribution
          - distribution_branches: Root branches used to access distribution
          - distribution_var_names: R variable names

    Returns:
        None: Run() stores a plot of the reconstructed spectrum
    """

    def __init__(self, name, *args, **kwargs):
        self.__name = name
        return

    def Configure(self, params):
        self.output_dir = read_param(params, 'output_dir', 'required')
        self.output_path_prefix = read_param(params, 'output_path_prefix', 'required')

        self.plot_data = read_param(params, 'plot_data', True)
        self.individual_spectra = read_param(params, 'make_individual_spectra', True)
        self.stacked_spectra = read_param(params, 'make_stacked_spectra', True)
        self.unstacked_spectra = read_param(params, 'make_unstacked_spectra', True)
        self.reconstruction_plot = read_param(params, 'make_reconstruction_plot', True)
        self.residual_plot = read_param(params, 'make_residual_plot', True)
        self.make_data_plot = read_param(params, 'make_data_plot', True)

        self.binning_file = read_param(params, 'binning_file', None)
        self.binning_file_format = read_param(params, 'binning_file_format', 'text')
        self.binning_file_variable = read_param(params, 'binning_file_variable', ':')
        self.n_bins = read_param(params, 'n_bins', 50)
        self.divide_by_bin_width = read_param(params, 'divide_by_bin_width', True)

        self.xmin = read_param(params, 'xmin', None)
        self.xmax = read_param(params, 'xmax', None)
        self.xlabel = read_param(params, 'xlabel', 'x')
        self.xlog = read_param(params, 'xlog', False)
        self.ybounds = read_param(params, 'ybounds', "auto")
        self.ylabel = read_param(params, 'ylabel', "y")
        self.ylog = read_param(params, 'ylog', False)
        self.title_prefix = read_param(params, 'title_prefix', "")

        if(self.plot_data or self.make_data_plot):
            self.data_path = \
                read_param(params, 'data_path', 'required')
            self.data_format = \
                read_param(params, 'data_format', 'required').split()
            self.data_variables = {}
            if "text" in self.data_format:
                self.data_variables["columns"] = \
                    read_param(params, 'data_columns', 'required')
            if "root" in self.data_format:
                self.data_variables["tree"] = \
                    read_param(params, 'data_tree', 'required')
                self.data_variables["branches"] = \
                    read_param(params, 'data_branches', 'required')
            if "R" in self.data_format:
                self.data_variables["variable_names"] = \
                    read_param(params, 'data_var_names', 'required')
            if "python" in self.data_format:
                self.data_variables["path"] = self.data_path
                self.data_variables["module"] = \
                    read_param(params, 'data_module', 'required')
                self.data_variables["method_name"] = \
                    read_param(params, 'data_function', 'required')
                self.data_variables["method_options"] = \
                    read_param(params, 'data_function_options', 'required')

        if(self.individual_spectra or self.stacked_spectra or
           self.unstacked_spectra):
            self.reconstructed_param_dicts = read_param(params, 'parameters', 'required')
            for p in self.reconstructed_param_dicts:
                p["shape_variables"] = {}
                if "text" in p["shape_format"]:
                    p["shape_variables"]["columns"] = \
                        read_param(p, 'shape_columns', 'required')
                if "root" in p["shape_format"]:
                    p["shape_variables"]["tree"] = \
                        read_param(p, 'shape_tree', 'required')
                    p["shape_variables"]["branches"] = \
                        read_param(p, 'shape_branches', 'required')
                if "R" in p["shape_format"]:
                    p["shape_variables"]["variable_names"] = \
                        read_param(p, 'shape_var_names', 'required')
                if "python" in p["shape_format"]:
                    p["shape_variables"]["path"] = p["shape_path"]
                    p["shape_variables"]["module"] = \
                        read_param(p, 'shape_module', 'required')
                    p["shape_variables"]["method_name"] = \
                        read_param(p, 'shape_function', 'required')
                    p["shape_variables"]["method_options"] = \
                        read_param(p, 'shape_function_options', 'required')
                p["distribution_variables"] = {}
                if "text" in p["distribution_format"]:
                    p["distribution_variables"]["columns"] = \
                        read_param(p, 'distribution_columns', 'required')
                if "root" in p["distribution_format"]:
                    p["distribution_variables"]["tree"] = \
                        read_param(p, 'distribution_tree', 'required')
                    p["distribution_variables"]["branches"] = \
                        read_param(p, 'distribution_branches', 'required')
                if "R" in p["distribution_format"]:
                    p["distribution_variables"]["variable_names"] = \
                        read_param(p, 'distribution_var_names', 'required')
                if "python" in p["distribution_format"]:
                    p["distribution_variables"]["path"] = p["distribution_path"]
                    p["distribution_variables"]["module"] = \
                        read_param(p, 'distribution_module', 'required')
                    p["distribution_variables"]["method_name"] = \
                        read_param(p, 'distribution_function', 'required')
                    p["distribution_variables"]["method_options"] = \
                        read_param(p, 'distribution_function_options', 'required')

        # Set up binnning
        try:
            self.binning = get_variable_from_file(self.binning_file,
                                                  self.binning_file_format,
                                                  self.binning_file_variable)
        except Exception as e:
            try:
                self.binning = np.linspace(self.xmin, self.xmax, self.n_bins+1)
            except Exception as e:
                print("Binning could not be created. Received error:")
                print(e)
                print("Exiting")
        self.bin_widths = []
        self.bin_centers = []
        for i in range(len(self.binning)-1):
            self.bin_widths.append(self.binning[i+1]-self.binning[i])
            self.bin_centers.append(self.binning[i] +
                                    (self.binning[i+1]-self.binning[i])/2.0)
        self.bin_widths = np.array(self.bin_widths, dtype='float64')
        self.bin_centers = np.array(self.bin_centers, dtype='float64')

        # Get data shape
        if(self.plot_data or self.make_data_plot):
            self.data_shape = self._get_histogram(self.data_path, self.data_format,
                                                  self.data_variables)[0]

        # Get parameter shapes and distributions
        if(self.individual_spectra or self.stacked_spectra or
           self.unstacked_spectra):
            for p in self.reconstructed_param_dicts:
                p["shape"] = self._get_histogram(p["shape_path"],
                                                 p["shape_format"],
                                                 p["shape_variables"])[0]
                if "function" in p["distribution_format"]:
                    logger.error("Format %s invalid, function not implemented for distributions"
                                 %p["distribution_format"])
                    logger.error("Distribution mean and average will not be accurate")
                temp_distribution = get_histo_shape_from_file([-float("inf"), float("inf")],
                                                              p["distribution_path"],
                                                              p["distribution_format"],
                                                              p["distribution_variables"])
                p["distribution_average"] = temp_distribution[1]
                p["distribution_sigma"] = temp_distribution[2]
        return

    def _get_histogram(self, path, file_format, variables):
        """Get a histogram from file, with the correct y axis

        If divide_by_bin_width is True, then we want histograms where
        the y axis is the number of counts per unit x. If False, then
        we want total number of counts as the y axis. This method
        uses the given path, file_format, and variables to get the
        histogram, then modifies it to have the correct y axis.

        Args:
            See morpho.file_reader.get_histo_shape_from_file docstring

        Returns:
            list of length (len(binning)-1), specifying either the
            counts per unit x for each bin, or the total number of
            counts in each bin.
        """
        (histo, avg, sigma) = get_histo_shape_from_file(self.binning, path,
                                                        file_format, variables)
        if not len(self.bin_widths)==len(histo):
            logger.warn("Histogram obtained from %s has invalid length"%path)
            logger.warn("Binning: %s"%self.binning)
            logger.warn("Bin widths: %s"%self.bin_width)
            logger.warn("Histogram: %s"%histo)
            logger.warn("Returning histogram")
            return histo

        if self.divide_by_bin_width:
            if "counts" in file_format or "values" in file_format:
                histo = histo/self.bin_widths
        else:
            if "function" in file_format:
                histo = histo*self.bin_widths
        return (histo, avg, sigma)

    def Run(self):
        """Create the plots"""
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if(self.divide_by_bin_width):
            histo_plot_type = "histo_line"
        else:
            histo_plot_type = "histo_points"

        if self.make_data_plot:
            data_curve = (self.binning, self.data_shape, histo_plot_type, {})
            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "data.png"
            plot_args = {"alpha":1.0}
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            plot_curves([data_curve], output_path, "matplotlib",
                        xlabel=self.xlabel, ylabel=self.ylabel,
                        title="Data",
                        xlog=self.xlog, ylog=self.ylog,
                        **plot_args)

        if self.individual_spectra:
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            for i,p in enumerate(self.reconstructed_param_dicts):
                curves = []
                if self.plot_data:
                    curves.append((self.binning, self.data_shape,
                                   histo_plot_type, {"label":"Data"}))
                curves.append((self.binning,
                               p["distribution_average"]*p["shape"],
                               histo_plot_type, {"label":p["name"]}))
                output_path = self.output_dir + "/" + \
                              self.output_path_prefix + "%s.png"%p["name"]
                plot_args = {}
                if not self.ybounds=="auto":
                    plot_args['ybounds'] = self.ybounds
                plot_curves(curves, output_path, plotter="matplotlib",
                            xlabel=self.xlabel, ylabel=self.ylabel,
                            title="Spectrum %s"%p["name"],
                            xlog=self.xlog, ylog=self.ylog,
                            **plot_args)

        if self.stacked_spectra:
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            curves = []
            total = np.zeros(len(self.binning)-1)
            for i,p in enumerate(self.reconstructed_param_dicts):
                total += p["distribution_average"]*p["shape"]
                curves = [(self.binning, copy.deepcopy(total),
                          "histo_shaded",
                           {"label":p["name"]})] + \
                          curves
            if self.plot_data:
                curves.append((self.binning, self.data_shape,
                               "histo_line",
                               {"label":"Data", "color":"black", "linewidth":2, "alpha":1}))
            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "stacked_spectra.png"
            plot_args = {"alpha":0.5}
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            plot_curves(curves, output_path, plotter="matplotlib",
                        xlabel=self.xlabel, ylabel=self.ylabel,
                        title="Stacked Spectra",
                        xlog=self.xlog, ylog=self.ylog,
                        **plot_args)

        if self.unstacked_spectra:
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            curves = []
            for i,p in enumerate(self.reconstructed_param_dicts):
                curves = [(self.binning,
                           p["distribution_average"]*p["shape"],
                          "histo_shaded",
                           {"label":p["name"]})] + \
                          curves
            if self.plot_data:
                curves.append((self.binning, self.data_shape,
                               "histo_line",
                               {"label":"Data", "color":"black", "linewidth":2, "alpha":1}))
            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "unstacked_spectra.png"
            plot_args = {"alpha":0.5}
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            plot_curves(curves, output_path, plotter="matplotlib",
                        xlabel=self.xlabel, ylabel=self.ylabel,
                        title="Unstacked Spectra",
                        xlog=self.xlog, ylog=self.ylog,
                        **plot_args)

        if self.reconstruction_plot:
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            curves = []
            total = np.zeros(len(self.binning)-1)
            for i,p in enumerate(self.reconstructed_param_dicts):
                total += p["distribution_average"]*p["shape"]
            curves.append((self.binning, total,
                           "histo_line",
                           {"label":"Reconstructed Spectrum"}))
            if self.plot_data:
                curves.append((self.binning, self.data_shape,
                               "histo_line",
                               {"label":"Data"}))
            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "reconstruction.png"
            plot_args = {"alpha":0.5}
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            plot_curves(curves, output_path, plotter="matplotlib",
                        xlabel=self.xlabel, ylabel=self.ylabel,
                        title="Reconstructed Spectrum",
                        xlog=self.xlog, ylog=self.ylog,
                        **plot_args)

        if self.residual_plot:
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            curves = []
            total = np.zeros(len(self.binning)-1)
            for i,p in enumerate(self.reconstructed_param_dicts):
                total += p["distribution_average"]*p["shape"]
            curves.append((self.binning, total,
                           "histo_line",
                           {"label":"Reconstructed Spectrum"}))
            if self.plot_data:
                curves.append((self.binning, self.data_shape,
                               "histo_line",
                               {"label":"Data"}))
            curves.append((self.binning,(self.data_shape-total),
                           "histo_line", {"label":"Residual (Exp-Stan)"}))
            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "reconstruction.png"
            plot_args = {"alpha":0.5}
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            plot_curves(curves, output_path, plotter="matplotlib",
                        xlabel=self.xlabel, ylabel=self.ylabel,
                        title="Reconstructed Spectrum",
                        xlog=self.xlog, ylog=self.ylog,
                        **plot_args)
        return


def reconstructed_spectrum(param_dict):
    """Create a spectrum reconstructed from stan parameters

    Creates a ReconstructSpectrumProcessor object, configures with
    param_dict, and then runs

    Args:
        param_dict: dictionary used to configure the
            ReconstructSpectrumProcessor object

    Returns:
        None: Stores a plot of the reconstructed spectrum
    """
    proc = ReconstructSpectrumProcessor("reconstructed_spectrum_plotter")
    proc.Configure(param_dict)
    proc.Run()
    return
