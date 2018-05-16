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
import matplotlib.pyplot as plt

try:
    import ROOT
except ImportError as e:
    logger.info("ROOT could not be imported, error: %s"%e)
    logger.info("Continuing without ROOT")

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
        individual_param_output_dir: String specifying output directory
            for plots of individual parameters. A separate directory
            can be specified because there could be a lot of these plots.

        output_path_prefix: String specifying output path prefix
        output_format: Format of output plots (Default="png")

        store_param_dists: Whether the results of a gaussian
            fit of each parameter distribution should be stored. (Default=True)
        store_param_dists_dir: Directory where the results of the
            gaussian fit should be stored (Default=output_dir)

        make_individual_spectra: Boolean specifying whether plots
            should be saved with each spectrum plotted separately
            (Default=True)
        make_stacked_spectra: Boolean specifying whether a plot should
            be saved with all spectra on the same plot (Default=True)
        make_unstacked_spectra: Boolean specifying whether a plot should
            be saved with all spectra on the same plot, stacked such that
            the summed spectra should approximately add to the data.
            (Default=True)
        make_reconstruction_plot: Whether a plot of the reconstructed
            spectrum should be stored (Default=True)
        make_diff_plot: Whether a plot of the data, reconstructed
            spectrum, and the difference should be stored (Default=True)
        make_residual_pull_plot: Whether a plot of normalized residuals
            and the pull histogram should be stored. (Default=True)
        make_data_model_ratio_plot: Whether a plot of the data to model
            ratio should be stored. (Default=True)
        make_chi2_vs_dof_plot: Whether a plot of the chi2 value vs the
            degrees of freedom should be stored. (Default=True)
        make_data_plot: Whether a plot should be made containing only
            the data. (Default=True)
        make_param_dist_plots: Whether a parameter distribution plot
            should be saved for each parameter. (Default=True)

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
        self.individual_param_output_dir = read_param(params,
                                                      'individual_param_output_dir',
                                                      self.output_dir)


        self.output_path_prefix = read_param(params, 'output_path_prefix', 'required')
        self.output_format = read_param(params, 'output_format', 'png')

        self.store_param_dists = \
            read_param(params, 'store_param_dists', True)
        self.store_param_dists_dir = read_param(params, 'param_dists_dir',
                                                  self.output_dir+"/param_dists")

        self.individual_spectra = read_param(params, 'make_individual_spectra', True)
        self.stacked_spectra = read_param(params, 'make_stacked_spectra', True)
        self.unstacked_spectra = read_param(params, 'make_unstacked_spectra', True)
        self.reconstruction_plot = read_param(params, 'make_reconstruction_plot', True)
        self.diff_plot = read_param(params, 'make_diff_plot', True)
        self.residual_pull_plot = read_param(params, 'make_residual_pull_plot', True)
        self.data_model_ratio_plot = read_param(params, 'make_data_model_ratio_plot', True)
        self.chi2_vs_dof_plot = read_param(params, 'make_chi2_vs_dof_plot', True)
        self.make_data_plot = read_param(params, 'make_data_plot', True)
        self.param_dist_plots = read_param(params, 'make_param_dist_plots', True)

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
                logger.error("Binning could not be created. Received error:")
                logger.error(e)
                logger.error("Exiting")
        self.bin_widths = []
        self.bin_centers = []
        for i in range(len(self.binning)-1):
            self.bin_widths.append(self.binning[i+1]-self.binning[i])
            self.bin_centers.append(self.binning[i] +
                                    (self.binning[i+1]-self.binning[i])/2.0)
        self.bin_widths = np.array(self.bin_widths, dtype='float64')
        self.bin_centers = np.array(self.bin_centers, dtype='float64')

        # Get data shape
        data_info = self._get_histogram(self.data_path, self.data_format,
                                        self.data_variables)
        self.data_shape = data_info[0]
        self.data_errors = data_info[3]

        # Get parameter shapes, distributions, and the total reconstruction
        self.tot_recon = np.zeros(len(self.binning)-1)
        self.tot_recon_errors = np.zeros(len(self.binning)-1) # sum errors in quadrature
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
            self.tot_recon += p["shape"]*p["distribution_average"]
            self.tot_recon_errors += (p["shape"]*p["distribution_sigma"])**2
        self.tot_recon_errors = np.sqrt(self.tot_recon_errors)

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
            4-tuple. First element is a list of length
            (len(binning)-1), specifying either the
            counts per unit x for each bin, or the total number of
            counts in each bin. Second is the average, third is sigma,
            fourth is the poisson error bar in each bin.
        """
        (histo, avg, sigma) = get_histo_shape_from_file(self.binning, path,
                                                        file_format, variables)
        if not len(self.bin_widths)==len(histo):
            logger.warn("Histogram obtained from %s has invalid length"%path)
            logger.warn("Binning: %s"%self.binning)
            logger.warn("Bin widths: %s"%self.bin_width)
            logger.warn("Histogram: %s"%histo)
            logger.warn("Returning histogram")
            return (histo, avg, sigma, histo*0.0)

        if self.divide_by_bin_width:
            if "counts" in file_format or "values" in file_format:
                histo = histo/self.bin_widths
            histo_pois = np.sqrt(histo*self.bin_widths)/self.bin_widths
        else:
            if "function" in file_format:
                histo = histo*self.bin_widths
            histo_pois = np.sqrt(histo)
        return (histo, avg, sigma, histo_pois)

    def Run(self):
        """Create the plots"""
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.exists(self.individual_param_output_dir):
            os.makedirs(self.individual_param_output_dir)
        if not os.path.exists(self.store_param_dists_dir):
            os.makedirs(self.store_param_dists_dir)

        if(self.divide_by_bin_width):
            histo_plot_type = "histo_line"
        else:
            histo_plot_type = "histo_points"

        if self.make_data_plot:
            data_curve = (self.binning, self.data_shape, "histo_error",
                          {"yerr":self.data_errors})
            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "data." + \
                          self.output_format
            plot_args = {"alpha":1.0}
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            plot_curves([data_curve], output_path, "matplotlib",
                        xlabel=self.xlabel, ylabel=self.ylabel,
                        title=self.title_prefix+"Data",
                        xlog=self.xlog, ylog=self.ylog,
                        **plot_args)

        if self.individual_spectra:
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            for i,p in enumerate(self.reconstructed_param_dicts):
                curves = []
                curves.append((self.binning, self.data_shape,
                               "histo_error",
                               {"label":"Data", "yerr":self.data_errors,
                                "color":"blue"}))
                curves.append((self.binning,
                               p["distribution_average"]*p["shape"],
                               "histo_error",
                               {"label":p["name"],
                                "yerr":p["distribution_sigma"]*p["shape"],
                                "color":"red"}))
                output_path = self.individual_param_output_dir + "/" + \
                              self.output_path_prefix + \
                              "%s.%s"%(p["name"],self.output_format)
                plot_args = {}
                if not self.ybounds=="auto":
                    plot_args['ybounds'] = self.ybounds
                plot_curves(curves, output_path, plotter="matplotlib",
                            xlabel=self.xlabel, ylabel=self.ylabel,
                            title="%sSpectrum %s"%(self.title_prefix,p["name"]),
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
                           {})] + \
                          curves
            curves.append((self.binning, self.data_shape,
                           "histo_line",
                           {"label":"Data", "color":"black", "linewidth":2, "alpha":1}))
            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "stacked_spectra." + \
                          self.output_format
            legend_path = self.output_dir + "/" + \
                          self.output_path_prefix + "stacked_spectra_legend." + \
                          self.output_format
            plot_args = {"alpha":0.5}
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            plot_curves(curves, output_path, plotter="matplotlib",
                        xlabel=self.xlabel, ylabel=self.ylabel,
                        title=self.title_prefix+"Stacked Spectra",
                        xlog=self.xlog, ylog=self.ylog,
                        legend_loc=0,
                        legend_output_path=legend_path,
                        **plot_args)

        if self.unstacked_spectra:
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            curves = []
            for i,p in enumerate(self.reconstructed_param_dicts):
                curves = [(self.binning,
                           p["distribution_average"]*p["shape"],
                          "histo_shaded",
                           {})] + \
                          curves
            curves.append((self.binning, self.data_shape,
                           "histo_line",
                           {"label":"Data", "color":"black", "linewidth":2, "alpha":1}))
            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "unstacked_spectra." + \
                          self.output_format
            legend_path = self.output_dir + "/" + \
                          self.output_path_prefix + "unstacked_spectra_legend." + \
                          self.output_format
            plot_args = {"alpha":0.5}
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            plot_curves(curves, output_path, plotter="matplotlib",
                        xlabel=self.xlabel, ylabel=self.ylabel,
                        title=self.title_prefix+"Unstacked Spectra",
                        xlog=self.xlog, ylog=self.ylog,
                        legend_loc=0,
                        legend_output_path=legend_path,
                        **plot_args)

        if self.reconstruction_plot:
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            curves = []
            curves.append((self.binning, self.data_shape,
                           "histo_error",
                           {"label":"Data",
                            "yerr":self.data_errors, "color":"blue"}))
            curves.append((self.binning, self.tot_recon,
                           "histo_error",
                           {"label":"Reconstruction",
                            "yerr":self.tot_recon_errors, "color":"red"}))
            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "reconstruction." + \
                          self.output_format
            plot_args = {"alpha":0.5}
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            plot_curves(curves, output_path, plotter="matplotlib",
                        xlabel=self.xlabel, ylabel=self.ylabel,
                        title=self.title_prefix+"Reconstruction",
                        xlog=self.xlog, ylog=self.ylog,
                        **plot_args)

        if self.diff_plot:
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            curves = []
            curves.append((self.binning, self.data_shape,
                           "histo_error",
                           {"label":"Data",
                            "yerr":self.data_errors, "color":"blue"}))
            curves.append((self.binning, self.tot_recon,
                           "histo_error",
                           {"label":"Reconstruction",
                            "yerr":self.tot_recon_errors, "color":"red"}))
            curves.append((self.binning,(self.data_shape-self.tot_recon),
                           "histo_error",
                           {"label":"Diff (Data-Recon)", "color":"black",
                            "yerr":np.sqrt(self.tot_recon_errors**2+
                                           self.data_errors**2)}))
            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "residual." + \
                          self.output_format
            plot_args = {"alpha":0.5}
            if not self.ybounds=="auto":
                plot_args['ybounds'] = self.ybounds
            plot_curves(curves, output_path, plotter="matplotlib",
                        xlabel=self.xlabel, ylabel=self.ylabel,
                        title=self.title_prefix+"Reconstruction Difference",
                        xlog=self.xlog, ylog=False,
                        **plot_args)

        if self.residual_pull_plot:
            pulls = []

            plt.subplot(1, 2, 1)
            plot_args = {}
            for i_bin in range(len(self.bin_centers)):
                pulls.append((self.data_shape[i_bin]-self.tot_recon[i_bin])/
                             np.sqrt(self.tot_recon_errors[i_bin]**2+self.data_errors[i_bin]**2))
            curves = [((self.bin_centers), np.array(pulls), "default",
                       {"marker":"^", "color":"black", "linestyle":"None"})]
            plot_curves(curves, None, plotter="matplotlib",
                        xlabel=self.xlabel, ylabel="N Sigma",
                        title=self.title_prefix+"(Data-Recon)/Sigma",
                        xlog=False, ylog=False, subplot=True, **plot_args)

            plt.subplot(1, 2, 2)
            curves = []
            (pulls_hist, pulls_hist_edges) = \
                np.histogram(pulls, bins=20, range=(-4, 4))
            curves.append((pulls_hist_edges, pulls_hist, "histo_error",
                           {"yerr":np.sqrt(pulls_hist), "color":"black"}))
            def norm_gauss(x):
                return np.exp(-x**2/2.)/np.sqrt(2.*np.pi)
            norm_gauss = np.vectorize(norm_gauss)
            x_pts = np.linspace(pulls_hist_edges[0], pulls_hist_edges[-1])
            y_pts = norm_gauss(x_pts)*len(pulls)*\
                    (pulls_hist_edges[1]-pulls_hist_edges[0])
            curves.append((x_pts, y_pts, "default",
                           {"color":"red"}))
            plot_curves(curves, None, plotter="matplotlib",
                        xlabel="N Sigma", ylabel="N Counts",
                        title=self.title_prefix+"Pull Distribution",
                        xlog=False, ylog=False, subplot=True, **plot_args)

            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "residuals_pulls." + \
                          self.output_format
            plt.savefig(output_path)
            plt.close()

        if self.data_model_ratio_plot:
            plot_args = {}
            ax1 = plt.subplot(2, 1, 1)
            curves = []
            curves.append((self.binning, self.data_shape,
                           "histo_error",
                           {"label":"Data", "color":"blue",
                            "yerr":self.data_errors}))
            curves.append((self.binning, self.tot_recon,
                           "histo_error",
                           {"label":"Reconstruction", "color":"red",
                            "yerr":self.tot_recon_errors}))
            plot_curves(curves, None, plotter="matplotlib",
                        xlabel="", ylabel=self.ylabel,
                        title=self.title_prefix+"Reconstruction",
                        xlog=self.xlog, ylog=self.ylog,
                        subplot=True,
                        **plot_args)

            ax2 = plt.subplot(2, 1, 2, sharex=ax1)
            curves = []
            ratios = self.data_shape/self.tot_recon
            ratio_errors = ratios*np.sqrt((self.data_errors/self.data_shape)**2+
                                          (self.tot_recon_errors/self.tot_recon)**2)
            curves.append((self.binning, 6.*ratio_errors, "histo_shaded",
                           {"bottom":1-3.*ratio_errors,"color":"red"}))
            curves.append((self.binning, 4.*ratio_errors, "histo_shaded",
                           {"bottom":1-2.*ratio_errors,"color":"yellow"}))
            curves.append((self.binning, 2.*ratio_errors, "histo_shaded",
                           {"bottom":1-ratio_errors,"color":"cyan"}))
            curves.append((self.bin_centers, ratios, "default",
                           {"marker":".", "color":"black", "linestyle":"None"}))
            plot_curves(curves, None, plotter="matplotlib",
                        xlabel=self.xlabel, ylabel="Data/Model Ratio",
                        title="",
                        xlog=self.xlog, ylog=False,
                        subplot=True,
                        **plot_args)

            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "data_model_ratio." + \
                          self.output_format
            plt.savefig(output_path)
            plt.close()
            pass

        if(self.chi2_vs_dof_plot):
            #plot_args['xbounds'] = (0., len(self.data_shape))
            #plot_args['ybounds'] = (0., len(self.data_shape))
            curves = []
            chi2_vals = [0.]
            chi2 = 0.
            for i_bin in range(len(self.data_shape)):
                diff2 = (self.data_shape[i_bin]-self.tot_recon[i_bin])**2
                sigma2 = (self.tot_recon_errors[i_bin]**2+self.data_errors[i_bin]**2)
                if(np.isclose(sigma2,0.)):
                    if(np.isclose(diff2,0.)):
                        chi2 += 1.
                    else:
                        logger.error("chi2 cannot be calculated for bin %i, sigma2==%.3e, diff2==%.3e"%
                                     (i_bin, sigma2, diff2))
                        chi2 += float("inf")
                else:
                    chi2 += diff2/sigma2
                chi2_vals.append(chi2)
            chi2_vals = np.array(chi2_vals)
            dof = np.array(range(len(self.data_shape)+1))
            curves.append((dof, dof, "default", {"color":"red"}))
            curves.append((dof, chi2_vals, "default",
                           {"marker":"^", "markersize":3, "color":"black",
                            "alpha":1, "linestyle":"None"}))
            output_path = self.output_dir + "/" + \
                          self.output_path_prefix + "chi2_vs_dof." + \
                          self.output_format
            plot_curves(curves, output_path, plotter="matplotlib",
                        xlabel=self.xlabel, ylabel="Chi-2",
                        title=self.title_prefix+"Chi-Squared vs d.o.f.",
                        xlog=False, ylog=False, **plot_args)

        if self.param_dist_plots or self.store_param_dists:
            if self.store_param_dists:
                text_output_file = self.store_param_dists_dir + \
                                   "/" + self.output_path_prefix + \
                                   "param_distribution_gauss_fits.txt"
                param_dist_file = open(text_output_file, 'w')
                p_name_col_width = 40
                param_dist_file.write("P Num   " +
                                      "Parameter Name".ljust(p_name_col_width) + "\t"
                                      "Dist Mean   \t" +
                                      "Dist Sigma  \t" +
                                      "Fit Mu      \t" +
                                      "Fit Sigma   \t" +
                                      "mu>3*sig\t" +
                                      "5% Quantile \t" +
                                      "95% Quantile\t" +
                                      "90% Quantile\t" +
                                      "\n")

            for i_param, p in enumerate(self.reconstructed_param_dicts):
                if not p["distribution_format"]=="root values":
                    logger.notice("Plotting parameter distributions only "+
                                  "available for file format 'root values'")
                    logger.warn("Skipping parameter %s"%p[name])
                    continue
                else:
                    myfile = ROOT.TFile(p["distribution_path"], "READ")
                    tree = myfile.Get(p["distribution_variables"]["tree"])
                    branch = p["distribution_variables"]["branches"][0]
                    cut = p["distribution_variables"]["cut"]
                    # Setting lb>ub will force automatic bininng
                    nbins = 50
                    param_hist = ROOT.TH1F("param_hist", "", nbins, 1., -1.)
                    curr_counts = tree.Draw(branch+">>param_hist", cut, "goff")
                    def gaus_fit_fcn(x, p):
                        arg = (x[0]-p[1])/p[2]
                        if(x[0]>0):
                            return p[0]*np.exp(-0.5*arg**2)/np.sqrt(2*np.pi*p[2]**2)
                        else:
                            return 0.
                    gaus_fitter = ROOT.TF1("gaus_fitter", gaus_fit_fcn, -100, 100, 3)
                    gaus_fitter.SetParameters(curr_counts/float(nbins),
                                              p["distribution_average"],
                                              p["distribution_sigma"])
                    param_hist.Fit(gaus_fitter, "Q0", "goff")
                    gaus_norm = gaus_fitter.GetParameter(0)
                    gaus_mean = gaus_fitter.GetParameter(1)
                    gaus_sigma = gaus_fitter.GetParameter(2)
                    if (gaus_mean-3.*gaus_sigma)>=0.:
                        gaus_fit_good = True
                    else:
                        gaus_fit_good = False
                    quantile_fracs = np.array([0.05, 0.95, 0.9])
                    # quantiles will be filled by GetQuantiles
                    quantiles = np.zeros(len(quantile_fracs))
                    param_hist.GetQuantiles(len(quantile_fracs),
                                            quantiles, quantile_fracs)

                    if self.param_dist_plots:
                        binning = []
                        bin_centers = []
                        bin_contents = []
                        if gaus_fit_good:
                            cl_90pct_lb = quantiles[0]
                            cl_90pct_ub = quantiles[1]
                        else:
                            cl_90pct_lb = 0.
                            cl_90pct_ub = quantiles[2]
                        for i_bin in range(1,nbins+1):
                            binning.append(param_hist.GetBinLowEdge(i_bin))
                            bin_centers.append(param_hist.GetBinCenter(i_bin))
                            bin_contents.append(param_hist.GetBinContent(i_bin))
                            if binning[-1]<=cl_90pct_lb:
                                shaded_bin_low = i_bin
                            elif binning[-1]<cl_90pct_ub:
                                shaded_bin_high = i_bin
                        binning.append(param_hist.GetBinLowEdge(nbins+1))
                        binning = np.array(binning)
                        bin_centers = np.array(bin_centers)
                        bin_contents = np.array(bin_contents)

                        curves = []
                        # Append 90% CL region, shaded
                        if gaus_fit_good:
                            temp_opts = {"color":"yellow", "alpha":0.7}
                        else:
                            temp_opts = {"color":"orange", "alpha":0.7}
                        curves.append((binning[shaded_bin_low:shaded_bin_high+1],
                                       bin_contents[shaded_bin_low:shaded_bin_high],
                                       "histo_shaded", temp_opts))

                        # Append histogram with poisson error bars
                        lb = float(param_hist.GetBinLowEdge(1))
                        ub = float(param_hist.GetBinLowEdge(nbins+1))
                        bin_width = (ub-lb)/float(nbins)
                        curves.append((binning, bin_contents, "histo_error",
                                       {"yerr":np.sqrt(bin_contents),
                                        "color":"black"}))

                        # Append gaussian fit, in blue or red depending on gaus_fit_good
                        if gaus_fit_good:
                            temp_opts_g = {"color":"blue"}
                        else:
                            temp_opts_g = {"color":"red"}
                        xpts = np.linspace(lb, ub, 200)
                        def curr_gaussian(x):
                            return gaus_fit_fcn([x], [gaus_norm, gaus_mean, gaus_sigma])
                        curr_gaussian = np.vectorize(curr_gaussian)
                        ypts = curr_gaussian(xpts)
                        curves.append((xpts, ypts, "default", temp_opts_g))

                        plot_output_file = self.store_param_dists_dir + \
                                           "/" + self.output_path_prefix + \
                                           "p%i_%s.%s"%(i_param, p["name"],
                                                         self.output_format)
                        plot_curves(curves, plot_output_file, plotter="matplotlib",
                                    xlabel=p["name"], ylabel="Distribution",
                                    title=self.title_prefix+p["name"]+" Distribution",
                                    xlog=False, ylog=False)

                    if self.store_param_dists:
                        param_dist_file.write("%i\t%s\t%.5e\t%.5e\t%.5e\t%.5e\t%s\t%.5e\t%.5e\t%.5e\n"%
                                              (i_param+1,
                                               p["name"][:p_name_col_width-2].ljust(p_name_col_width, '.'),
                                               p["distribution_average"],
                                               p["distribution_sigma"],
                                               gaus_mean,
                                               gaus_sigma,
                                               ("%s"%gaus_fit_good).ljust(8),
                                               quantiles[0],
                                               quantiles[2],
                                               quantiles[1]))
            param_dist_file.close()
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
