#======================================================
# list_plotter
#
# Author: J. Johnston
# Date: Mar. 1, 2018
#
# Make plots from lists
#=======================================================

"""Make plots from python lists

Functions:
  - mpl_plot_curves: Plot curves using matplotlib
  - plot_curves: Plot curves with various plotting modules
"""

import logging
logger = logging.getLogger(__name__)
import copy

import numpy as np

import matplotlib as mpl
mpl.rc('ytick', labelsize=8)
mpl.rc('xtick', labelsize=8)
import matplotlib.pyplot as plt

def mpl_plot_curves(curves, output_path,
                    xlabel="", ylabel="",title="",
                    xbounds=None, xlog=False,
                    ybounds=None, ylog=False,
                    colors=['black', 'blue', 'green', 'red', 'cyan',
                            'magenta', 'yellow', 'orange', 'purple'],
                    alpha=1.0, legend_size=None, legend_loc=0,
                    legend_output_path=None):
    """Plot curves using matplotlib

    Args:
        curves: A list of 4-tuples, with each4-tuple containing
            (x_pts, y_pts, curve_type, opts)
            where possible curve_type options are, "default",
            "histo_shaded", "histo_line", "histo_points", and
            "histo_error". Default simply uses
            plt.plot to plot the given points. "histo_shaded", etc
            require x_pts to have length one greater than
            y_pts, where x_pts define the bin edges, and y_pts specify
            the bin contents. "histo_points" also plots a point at the
            top and center of each bin, and "histo_error" also plots
            error bars. In that case, y_pts should be 2xN, where the first
            row gives the y points, and the second row gives the error bar.
            opts is a dictionary of plotting options,
            such as a legend label, and marker and line style. If "color"
            or "alpha" is defined in opts, this will override the given
            colors and alpha.
        output_path: Location to store the resulting plot
        xlabel, ylabel, title: Labels for the plot
        xbounds, ybounds: 2-tuples with (lower_bound, upper_bound)
        xlog, ylog: Whether the axis should be log scaled
        colors: Colors used to plot the curves
        alpha: Transparency of the curves
        legend_size: int giving legend size

    Returns:
        None: Stores a plot at the given path.
    """
    figData = plt.figure()

    for i,c in enumerate(curves):
        try:
            xpts = c[0]
            ypts = c[1]
            c_type = c[2]
            opts = c[3]
            if not "color" in opts:
                opts["color"] = colors[i%len(colors)]
            if not "alpha" in opts:
                opts["alpha"] = alpha

            if(c_type=="histo_shaded"):
                bar_widths = []
                for j in range(len(xpts)-1):
                    bar_widths.append(xpts[j+1]-xpts[j])
                if not "linewidth" in opts:
                    opts["linewidth"] = 0
                plt.bar(xpts[:-1], ypts,
                        width=bar_widths, **opts)
            elif(c_type=="histo_line"):
                histo_ypts = np.append(min(ypts), ypts)
                plt.plot(xpts, histo_ypts, ls='steps', **opts)
            elif(c_type=="histo_points"):
                bin_centers = []
                for i in range(len(xpts)-1):
                    bin_centers.append(xpts[i]+
                                       (xpts[i+1]-xpts[i])/2.0)
                histo_ypts = np.append(min(ypts), ypts)
                temp_opts = copy.deepcopy(opts)
                temp_opts["marker"] = "None"
                plt.plot(xpts, histo_ypts, ls='steps', **temp_opts)
                temp_opts = copy.deepcopy(opts)
                temp_opts["linestyle"] = "None"
                if not "marker" in temp_opts:
                    temp_opts["marker"] = '.'
                plt.plot(bin_centers, ypts, **temp_opts)
            elif(c_type=="histo_error"):
                bin_centers = []
                for i in range(len(xpts)-1):
                    bin_centers.append(xpts[i]+
                                       (xpts[i+1]-xpts[i])/2.0)
                histo_ypts = np.append(min(ypts), ypts)
                temp_opts = copy.deepcopy(opts)
                temp_opts["marker"] = "None"
                plt.plot(xpts, histo_ypts, ls='steps', **temp_opts)
                temp_opts = copy.deepcopy(opts)
                temp_opts["linestyle"] = "None"
                if not "marker" in temp_opts:
                    temp_opts["marker"] = '*'
                plt.plot(bin_centers, ypts, **temp_opts)
                print("Plotting a histo with error bars not yet implemented- Just plots with *")
            else:
                if not c_type=="default":
                    logger.info("Invalid curve type. Using plt.plot")
                plt.plot(xpts, ypts, **opts)
        except Exception as e:
            logger.warn("While plotting curve %i, got exception %s"%(i,e))

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if(xlog):
        plt.gca().set_xscale('log')
    if(ylog):
        plt.gca().set_yscale('log')

    legend_dict = {}
    if(not legend_size is None):
        legend_dict['prop'] = {'size':legend_size}
    if(legend_loc is None):
        pass
        # Need to update this to save the legend to legend_output_path
        '''ax = plt.gca()
        figData.savefig(output_path)
        figData.clf()
        handles, leg_labels = ax.get_legend_handles_labels()
        fig_legend = plt.figure(figsize=(2,2))
        axi = fig_legend.add_subplot(111)
        fig_legend.legend(handles, leg_labels, loc='center', scatterpoints = 1)
        axi.xaxis.set_visible(False)
        axi.yaxis.set_visible(False)
        if not legend_output_path is None:
            fig_legend.savefig(legend_output_path)
        else:
            fig_legend.savefig(output_path+"_legend.png")
        fig_legend.clf()'''
    else:
        legend_dict['loc'] = legend_loc
        plt.legend(**legend_dict)
    figData.savefig(output_path)
    figData.clf()
    plt.close()
    return

def plot_curves(curves, output_path,
                plotter="matplotlib",
                xlabel="", ylabel="",title="",
                xbounds=None, xlog=False,
                ybounds=None, ylog=False,
                colors=['black', 'blue', 'green', 'red', 'cyan',
                        'magenta', 'yellow', 'orange', 'purple'],
                alpha=1.0, legend_size=None,
                legend_loc=0, legend_output_path=None):
    """Plot a histogram using matplotlib

    Args:
        curves: A list of 4-tuples, with each4-tuple containing
            (x_pts, y_pts, curve_type, opts)
            where possible curve_type options are, "default",
            "histo_shaded", "histo_line", "histo_points", and
            "histo_error". Default simply uses
            plt.plot to plot the given points. "histo_shaded", etc
            require x_pts to have length one greater than
            y_pts, where x_pts define the bin edges, and y_pts specify
            the bin contents. "histo_points" also plots a point at the
            top and center of each bin, and "histo_error" also plots
            error bars. In that case, y_pts should be 2xN, where the first
            row gives the y points, and the second row gives the error bar.
            opts is a dictionary of plotting options,
            such as a legend label, and marker and line style. If "color"
            or "alpha" is defined in opts, this will override the given
            colors and alpha.
        output_path: Location to store the resulting plot
        plotter: Specify plotter to use. Options are "matplotlib" and
            "root" ("root" is not yet implemented).
        xlabel, ylabel, title: Labels for the plot
        xbounds, ybounds: 2-tuples with (lower_bound, upper_bound)
        xlog, ylog: Whether the axis should be log scaled
        colors: Colors used to plot the curves
        alpha: Transparency of the curves
        legend_size: int giving legend size
    """
    if plotter=="root":
        logger.warn("root plotting not yet implemented. Using matplotlib.")
    else:
        # Use matplotlib
        pass
    mpl_plot_curves(curves, output_path, xlabel, ylabel, title,
                    xbounds, xlog, ybounds, ylog,
                    colors, alpha, legend_size,
                    legend_loc, legend_output_path)
    return
