#!/usr/bin/evn python
#
# Class to perform a frequentist sensitivity analysis
# ----------------------------------------------------------------------
# Author: J. P. Johnston <jpj13@mit.edu>
#
# ----------------------------------------------------------------------

import numpy as np
import scipy as sp

import ROOT
import root_tools as rt

from array import array
import os

try:
    from pystan.misc import stan_rdump
except ImportError, e:
    pass

class InvalidBoundsError(Exception):
    """Throw this when bounds are invalid"""

class DimensionMismatchError(Exception):
    """Throw this when there is not a shape in a Parameter object matching
    the current Dimension in an Experiment"""

class MinimizationFailedError(Exception):
    """Throw this when minimization is not successful"""
    
class Dimension:
    """Data structure describing one dimension of a sensitivity analysis"""

    def set_bounds(self,lb,ub):
        self._lb = lb
        if(ub>lb):
            self._ub = ub
        else:
            self._ub = lb
            raise InvalidBoundsError("upper bound < lower bound, setting ub=lb")
        return

    def get_bounds(self):
        """Returns length 2 list (lb,ub)"""
        return (self._lb,self._ub)

    def __init__(self, name, lb, ub, nbins=10,
                 multiply_rate=False, convolution_fcn=None):
        """name: name of the dimension, eg time or energy
        lb, ub: Lower and upper bounds of the dimension
        nbins: Number of bins to use for a binned analysis
        multiply_rate: Boolean specifying whether rates should be multiplied by interval
                       length along this dimension. That is, if rates are in counts per second,
                       then multiply_rate should be true for the time dimensions.
        convolution_fcn: Function to convolve with shapes in order to account for resolution"""
        self.name = name
        self.nbins = nbins
        self.multiply_rate = multiply_rate
        self.convolution_fcn = convolution_fcn
        self.set_bounds(lb,ub)
        return
    
    def __repr__(self):
        str = "Name=%s, " % self.name
        str += "Bounds=[%.2e,%.2e], " % self.get_bounds()
        str += "NumBins=%i, " % self.nbins
        str += "Multiply_Rate=%s" % self.multiply_rate
        return str

    def get_bin_width(self):
        return (self._ub-self._lb)/float(self.nbins)

    
class Shape:
    """Data structure to store a shape for a sensitivity analysis"""

    def __init__(self, name, fcn, fcn_args=list(), 
                 dimension_name="", samples_per_bin=100):
        """
        fcn: Can either be a python function or a 2d numpy array.
             If fcn is a python fcn, then fcn(x,fcn_args) should take a numpy array
             x and return a numpy array of values y
             If fcn is an array then x values should be in fcn[0,:] and y values 
             should be in fcn[1,:]. Note that two 1d arrays can be
             combined in this way with fcn = np.vstack((xvals,yvals)).
        fcn_args: Arguments to pass. Use only if fcn is a python function.
        name: String used to identify this shape when printing messages.
        dimension_name: Name associated with the Dimension object corresponding to 
                        this shape. Used by Experiment to match shapes with dimensions
        samples_per_bin: The number of equally spaced points to sample in order to get 
                         the average y value in each bin. Only used in the case of a 
                         binned analysis where fcn is python fcn.
        """
        self.fcn = fcn
        self.fcn_args = fcn_args
        self.name = name
        self.dimension_name = dimension_name
        self.samples_per_bin = samples_per_bin
        self._shape_saved_loc = None # Last place the shape was saved
        return

    def __repr__(self):
        str = "Shape %s: " % (self.name)
        str += "fcn=%s, " % self.fcn
        str += "fcn_args=%s, " % self.fcn_args
        str += "dimension:%s, " % self.dimension_name
        str += "samples_per_bin=%i, " % self.samples_per_bin
        str += "last saved in %s" % self._shape_saved_loc
        return str

    def __str__(self):
        str = "%s: " % (self.name)
        str += "dimension=%s, " % self.dimension_name
        if(callable(self.fcn)):
            str += "fcn_type=py, "
            str += "samps_per_bin=%i, " % self.samples_per_bin
        else:
            str += "fcn_type=arrays, "
        str += "last saved in %s" % self._shape_saved_loc
        return str

    def generate_shape(self,lb, ub, renormalize=False, binned=True,
                       nbins=10, convolution_fcn=None):
        """Generate a shape with bounds lb and ub, optionally renormalized to 1
        If binned=True, then a 2D numpy array will be returned,
          with x values in shape[0,:] and y values in shape[1,:]. The x values
          are the value in the center of each bin, and the y values are the 
          average value of fcn in each bin.
        Otherwise, a python function of 1 variable will be returned.
        convolution_fcn must be a python function. If one is provided, then
          fcn will be convolved with convolution_fcn before generating the shape"""
        if(not convolution_fcn is None):
            # Figure out how to do a convolution
            print("Convolution not yet implemented")
            temp_fcn = self.fcn
        else:
            temp_fcn = self.fcn

        if(binned):
            if(callable(temp_fcn)):
                x_vals = np.linspace(lb,ub,nbins*self.samples_per_bin,False)
                if(len(self.fcn_args)>0):
                    y_vals = temp_fcn(x_vals,self.fcn_args)
                else:
                    y_vals = temp_fcn(x_vals)
            else:
                x_vals = temp_fcn[0,:]
                y_vals = temp_fcn[1,:]
            # Get the average value in each bin
            h = ROOT.TH1F("h","",nbins,lb,ub)
            hw = ROOT.TH1F("hw","",nbins,lb,ub)
            rebinned_x = []
            rebinned_y = []
            for j in range(0,len(x_vals)):
                h.Fill(x_vals[j],y_vals[j])
                hw.Fill(x_vals[j],1)
            # Bin numering starts at 1, not 0. List numbering starts at 0
            for j in range(1,nbins+1):
                rebinned_x.append(h.GetBinCenter(j))
                rebinned_y.append(h.GetBinContent(j)/max(1,hw.GetBinContent(j)))
            if(renormalize):
                rebinned_y = np.asarray(rebinned_y)/float(sum(rebinned_y))
            return np.vstack((rebinned_x,rebinned_y))
        else:
            print("Unbinned shape generation not yet implemented")
            print("\tReturning None")
            return None

    def save_shape_root(self,binned_shape, path_root, x_br_name="x",
                        y_br_name="y", tree_name=None, recreate_file=True):
        """Save the points for a binned shape in a root file
        binned_shape: Output of generate_shape where binned=True"""
        if(callable(binned_shape)):
            print("Shape is not binned and cannot be saved in a root file")
            return

        if(path_root[-5:]!=".root"):
            path_root += ".root"
        if(tree_name is None):
            tree_name=self.name.replace(" ","_")
        if(recreate_file):
            myfile = ROOT.TFile(path_root,"RECREATE")
        else:
            myfile = ROOT.TFile(path_root,"UPDATE")
        # Create the TTree
        tmp_x = array('f',[ 0 ])
        tmp_y = array('f',[ 0 ])
        shape_tree = ROOT.TTree(tree_name,tree_name)
        shape_tree.Branch(x_br_name,tmp_x,x_br_name+'/F')
        shape_tree.Branch(y_br_name,tmp_y,y_br_name+'/F')
        for j in range(0,len(binned_shape[0,:])):
            tmp_x[0] = binned_shape[0,j]
            tmp_y[0] = binned_shape[1,j]
            shape_tree.Fill()
        shape_tree.Write()
        myfile.Close()

        self._shape_saved_loc = "%s, tree=%s, x_br=%s, y_br=%s" % (path_root, tree_name, x_br_name, y_br_name)
        return

    def save_shape_pdf(self, shape, path_pdf,
                         x_label="x", y_label="y", title=None):
        """Plot the shape and save in a pdf file
        shape: output of generate_shape"""
        if(path_pdf[-4:]!=".pdf"):
            path_pdf += ".pdf"
        if(title is None):
            title = self.name

        if(shape is callable(shape)):
            print("Unbinned stuff needs to be implemented")
        else:
            temp_gr = ROOT.TGraph(len(shape[0,:]),np.array(shape[0,:]),np.array(shape[1,:]))
            rt.root_make_plot([temp_gr],path_pdf,["APL"],title,x_label,y_label,marker_styles=[3])
        self._shape_saved_loc = path_pdf
        return

class Parameter:
    """Data structure describing one parameter of a sensitivity analysis"""

    def __init__(self, name, shapes=list(), renormalize=False,
                 nuisance=False, gauss_constraint=False,
                 gauss_mean=None, gauss_sigma=None,
                 magnitude=0.0, magnitude_H0=None, lb=None, ub=None):
        """Initialize a parameter object with the following data:
        name: Identifier for the parameter
        shapes: list of Shape objects. Shapes can also be added with add_shape()
                In a sensitivity analysis, shapes should be added in the same order
                as the corresponding dimensions.
        renormalize: Whether shapes for this parameter should be renormalized to 1.
                     For example, if the parameter is a number of counts or a rate
                     then it should be renormalized to 1. But if the parameter is
                     a physical constant then it should not be renormalized.
        nuisance: Boolean defining whether this is a nuisance parameter
        gauss_constraint: Boolean defining whether this parameter is constrained
                          by external data. The data is assumed to be gaussian.
        gauss_mean: Mean for the gaussian contraint
        gauss_sigma: Sigma for the gaussian constraint
        magnitude: Current magnitude of this parameter
        lb, ub: Give bounds on this parameter when generating a value using 
                the gaussian contraint (ie constraint for H1 hypothesis)
        magnitude_H0: magnitude of this parameter under the H0 constraint. Set to
                      None if this parameter is not constrained under H0
        """
        self.name = name
        self.shapes = shapes
        self.renormalize = renormalize
        self.nuisance = nuisance
        self.gauss_constraint = gauss_constraint
        self.gauss_mean = gauss_mean
        self.gauss_sigma = gauss_sigma
        self.magnitude = magnitude
        self.magnitude_H0 = magnitude_H0
        self.lb = lb
        self.ub = ub
        self.extracted_tree = None
        return

    def __str__(self):
        str = "%s: " % self.name
        str += "Nuisance=%s, " % self.nuisance
        if(self.gauss_constraint):
            str += "Gauss Constrained, Mean=%.3e, Sigma=%.3e" % (self.gauss_mean,self.gauss_sigma)
            if(not self.lb is None):
                str += ", lb=%.2e" % self.lb
            if(not self.ub is None):
                str += ", ub=%.2e" % self.ub
            str += "; "
        else:
            str += "Not Gauss Constrained; "
        if(not self.magnitude_H0 is None):
            str += "H0_mag=%.2e; " % self.magnitude_H0
        str += "Shapes- "
        for s in self.shapes:
            str += "%s; " % s
        str += "Magnitude=%.3e" % self.magnitude
        return str

    def add_shape(self, shp):
        """Shapes should be added in the same order that Dimensions were added
        The array of shapes can be accessed with the variable shapes"""
        self.shapes.append(shp)
        return

    def generate_magnitude(self, store=True):
        """If the parameter is constrained by a gaussian, randomly generate a magnitude
        with that gaussian and return it. If store=True, then the mean of that gaussian
        in self.magnitude before returning.
        Otherwise return the value stored in self.magnitude
        Thus the value of self.magnitude for each parameter will be the true
        values , while the values in the returned array will be the gaussian distributed
        values."""
        if(self.gauss_constraint):
            temp_mag = np.random.normal(self.gauss_mean,self.gauss_sigma)
            if(store):
                self.magnitude = self.gauss_mean
        else:
            temp_mag = self.magnitude
        return temp_mag

    def log_likelihood_gauss_factor(self, magnitude=None):
        """If a parameter is constrained by a gaussian, then the total likelihood will
        be multiplied by the factor gauss(curr_val|constraint_mean,constraint_sigma).
        if self.gauss_constraint is true, then likeihood_gauss_factor returns that factor
        for the given magnitude or, if magnitude=None, the object's magnitude.
        Otherwise it returns 1.0"""
        if(self.gauss_constraint):
            if(magnitude is None):
                x = self.magnitude
            else:
                x = magnitude
            mu = self.gauss_mean
            sig = self.gauss_sigma
            return sp.stats.norm.logpdf(x,self.gauss_mean,self.gauss_sigma)
        else:
            return 1.0

    def extracted_tree_initialize(self):
        """Initialize the histogram that contains parameters.
        A root file must be opened before initializing, and it must remain open
        as long as you want to use histogram of extracted values"""
        # EDIT THIS
        self.extracted_tree = ROOT.TTree()
        return


class Experiment:
    """Class to represent one counting experiment
    Data:
    
    dimensions: list of Dimension objects
    parameters: list of Parameter objects
    expt_settings: Settings specific to the experiment, for use by tot_rate_fcn
    tot_rate_fcn: Function to calculate total rate from lists. The inputs must
                  be (parameters, expt_settings, dimensions)
    data: Data for use in calculating the likelihood"""

    def reset(self, name, dimensions=list(), parameters=list(), binned_expt=True,
              expt_settings=None, tot_rate_fcn=None, data=None):
        """"Return Experiment object to initial state
        See __init__ documentation for details about inputs"""
        self.name = name
        self.dimensions = dimensions
        self.parameters = parameters
        self.binned_expt = binned_expt
        self.expt_settings = expt_settings
        self.tot_rate_fcn = tot_rate_fcn
        self.data = data
        self._data_saved_loc = None
        self._xy_1d_shapes = None
        self._nd_shapes = None
        return

    def __init__(self, name, dimensions=list(), parameters=list(), binned_expt=True,
                 expt_settings=None, tot_rate_fcn=None, data=None):
        """name: string to describe this experiment
        dimensions: list of dimension objects for this experiment.
        parameters: list of parameter objects that can be used in this experiment.
                    Each Parameter object must contain a shape for every Dimension
                    in the Experiment object's dimensions array.
        binned_expt: Whether data and the likelihood should be binned or unbinned.
        expt_settings: Array of settings for experiment, to be used by tot_rate_fcn.
        tot_rate_fcn: Function to calculate total rate from lists. The inputs must
                  be (parameters, expt_settings, dimensions)
        data:  Data used to calculate the likelihood. Can later be set with generate_fake_data()."""
        self.reset(name, dimensions, parameters, binned_expt,
                   expt_settings, tot_rate_fcn, data)
        return

    def __str__(self):
        str = "Experiment %s\n" % self.name
        str += "  Dimensions:\n"
        for d in self.dimensions:
            str += "    %s\n" % d
        str += "  Parameters:\n"
        for p in self.parameters:
            str += "    %s\n" % p
        str += "  Data last saved in: %s" % self._data_saved_loc
        return str


    def get_dimension(self,name):
        """Returns the first dimension object with the given name"""
        for d in self.dimensions:
            if d.name==name:
                return d
                break
        else:
            print("Dimension \"%s\" not found by Experiment %s get_dimension()" % (name,self.name))
        return

    def remove_dimension(self,name):
        """Removes and returns the first dimension with the given name"""
        for i,d in enumerate(self.dimensions):
            if d.name==name:
                return self.dimensions.pop(i)
        else:
            print("Dimension \"%s\" not found by Experiment %s remove_dimension()" % (name,self.name))
            return None

    def print_dimensions(self):
        print("Dimensions stored in Experiment %s" % self.name)
        for d in self.dimensions:
            print("\t%s" % d)

    def get_parameter(self, name):
        """Get a Parameter object by name, so it can then be edited"""
        for p in self.parameters:
            if p.name==name:
                return p
                break
        else:
            print("Parameter \"%s\" not found by Experiment %s get_parameter()"%(name,self.name))  
        return

    def remove_parameter(self, name):
        """Remove a the first Parameter object with the given name and return it"""
        for i,p in enumerate(self.parameters):
            if p.name==name:
                return self.parameters.pop(i)
                break
        else:
            print("Parameter \"%s\" not found by Experiment %s remove_parameter()"%(name,self.name))
            return None

    def print_parameters(self):
        """Print information about all parameters, including names used for access with get_parameter()"""
        str = "Parameters:\n"
        for p in self.parameters:
            str += "  %s\n" % p
        print(str)
        return


    def generate_all_shapes(self):
        """Generate all shapes for the given dimensions and parameters
        This only needs to be done once for a given set of dimensions and parameters
        (specifically, shapes do not need to be regenerated when the magnitude is 
        changed for a parameter)
        Each parameter must contain a shape with dimension_name equal to the name
        of one of the objects in dimensions. The bounds, number of bins, etc from
        that dimension will be used to generate the shape"""
        self._xy_1d_shapes = list()
        #self._y_1d_shapes = list()
        self._nd_shapes = list()
        binned = self.binned_expt
        for i,p in enumerate(self.parameters):
            self._xy_1d_shapes.append(list())
            for d in self.dimensions:
                for s in p.shapes:
                    if(d.name == s.dimension_name):
                        (lb,ub) = d.get_bounds()
                        renormalize = p.renormalize and not d.multiply_rate
                        nbins = d.nbins
                        convolution_fcn = d.convolution_fcn
                        sh = s.generate_shape(lb, ub, renormalize, binned,
                                                  nbins, convolution_fcn)
                        self._xy_1d_shapes[i].append(sh)
                        break
                else:
                    str = "Shape [%s] did not match any dimensions in " % s
                    str += "Experiment %s with dimensions: " % self.name
                    for d in self.dimensions:
                        str += "[%s], " % d
                    str = str[:-2]
                    raise DimensionMismatchError(str)
            # Create an n dimensional array or function from the shapes for this parameter
            if(binned):
                nd_arr = np.array([1])
                for j,shxy in enumerate(self._xy_1d_shapes[i]):
                    sh = shxy[1,:]
                    temp_arr = np.zeros(nd_arr.shape+sh.shape)
                    for index_0,val_0 in np.ndenumerate(nd_arr):
                        for index_1,val_1 in np.ndenumerate(sh):
                            temp_arr.itemset(index_0+index_1,val_0*val_1)
                    nd_arr = temp_arr
                self._nd_shapes.append(nd_arr[0])
            else:
                def temp_fcn(dim_vals_arr):
                    val = 1.0
                    for j,sh in self._xy_1d_shapes[i]:
                        val *= sh(dim_vals_arr[j])
                    return val
                temp_fcn = np.vectorize(temp_fcn)
                self._nd_shapes.append(temp_fcn)
        return

    def save_all_shapes_root(self, path_root, recreate=True):
        """Save all shapes for all parameters in a rootfile with the given filenam
        The tree for each shape will be the Shape name, the the x branch will use 
        the corresponding Dimension name and the y branch will use the Parameter name
        (although any spaces will be replaced by "_" in the tree and branch names)
        """
        if(path_root[-5:]!=".root"):
            path_root += ".root"
        for i,p in enumerate(self.parameters):
            for j,d in enumerate(self.dimensions):
                shape = self._xy_1d_shapes[i][j]
                if(i>0 or j>0):
                    recreate = False
                for s in p.shapes:
                    if(d.name == s.dimension_name):
                        x_br = d.name.replace(" ","_")
                        y_br = p.name.replace(" ","_")
                        tree = s.name.replace(" ","_")
                        s.save_shape_root(shape, path_root, x_br, y_br, tree, recreate)
        return

    def save_all_shapes_pdf(self, path_directory):
        """Save pdf plots of all shapes in the given directory
        The filename and title will be the Shape name. The y label will use the
        Parameter name, and the x label will use the corresponding Dimension name"""
        if not os.path.exists(path_directory):
            os.makedirs(path_directory)
        for i,p in enumerate(self.parameters):
            for j,d in enumerate(self.dimensions):
                shape = self._xy_1d_shapes[i][j]
                if(i==0 and j==0):
                    recreate = True
                else:
                    recreate = False
                for s in p.shapes:
                    if(d.name == s.dimension_name):
                        x_lab = d.name.replace(" ","_")
                        y_lab = p.name.replace(" ","_")
                        title = s.name.replace(" ","_")
                        path_pdf = path_directory + "/" + title + ".pdf"
                        s.save_shape_pdf(shape, path_pdf, x_lab, y_lab, title)
        return


    def calculate_expected_counts(self, param_mags=None):
        if(param_mags is None):
            param_mags = np.zeros(len(self.parameters))
            for i,p in enumerate(self.parameters):
                param_mags[i] = p.magnitude
        elif len(param_mags) != len(self.parameters):
            str = "param_mags {%s} and self.parameters {%s} have different dimensions" % (param_mags,self.parameters)
            raise DimensionMismatchError(str)
        exp_counts_arr = np.zeros(self._nd_shapes[0].shape)
        for (index,elt) in np.ndenumerate(exp_counts_arr):
            temp_mags = np.zeros(len(param_mags))
            for i,nd in enumerate(self._nd_shapes):
                temp_mags[i] = param_mags[i]*nd.item(index)
            exp_rate = self.tot_rate_fcn(temp_mags, self.expt_settings, self.dimensions)
            for d in self.dimensions:
                if d.multiply_rate:
                    exp_rate *= d.get_bin_width()
            exp_counts_arr.itemset(index,exp_rate)
        return exp_counts_arr

    def generate_fake_data(self):
        """Generate fake data using the current dimensions and parameters
        generate_all_shapes() must be called after storing all dimensions and parameters
        update the magnitude in each parameters that is not constrained before generating fake data
        Fake data is stored as an array in the variable data
          For a binned analysis, it will be an n dimensional array whose values signify counts
          For an unbinned analysis, it will be an array whose elements are positions where the
          events were generated (eg times, energies, or 2D arrays with time and energy)"""
        param_mags = np.zeros(len(self.parameters))
        for (i,p) in enumerate(self.parameters):
            param_mags[i] = p.generate_magnitude()
        if(self.binned_expt):
            exp_counts = self.calculate_expected_counts(param_mags)
            self.data = np.random.poisson(exp_counts)
        else:
            print("Unbinned analysis fake  data generation not yet implemented")
        return

    def _get_1d_data(self, dim_num):
        # All parameters will have the same x data, so we choose parameter 0
        x_data = self._xy_1d_shapes[0][dim_num][0,:]
        y_data = np.zeros(len(x_data))
        num_dims = len(self.dimensions)
        index = np.zeros(num_dims, dtype=np.int8)
        for i in range(0,len(y_data)):
            index[dim_num] = i
            y_data[i] = self.data.item(tuple(index))
        return np.vstack((x_data,y_data))
        

    def save_data_root(self, path_root, tree_name_prefix="data", br_name="counts",
                       recreate=True):
        """Save current binned data in a root file
        For multidimensional data, it saves one row for each dimension, with the
        other dimensions set to the 0th index
        "_dimension_name" is appended to tree_name_prefix
        The dimension name is used as the branch name for a branch with the x values
        """
        if(not self.binned_expt):
            print("Unbinned data cannot be saved to root file. Current experiment:\n %s" % self)
            return
        if(self.data is None):
            print("Data not set, and cannot be saved. Current experiment:\n %s" % self)
            return
        if(path_root[-5:]!=".root"):
            path_root += ".root"
        temp_sh = Shape("temp_shape",None)
        for i,d in enumerate(self.dimensions):
            if(i>0):
                recreate = False
            x_br = d.name
            tree_name = tree_name_prefix + "_" + d.name
            temp_sh.save_shape_root(self._get_1d_data(i), path_root,
                                    x_br, br_name, tree_name, recreate)
        self._data_saved_loc = path_root
        return

    def save_data_r(self, path_r, var_name="counts",
                    append_to_file=None, data_type=None):
        """Save current data in the r file format
        path_r: File path to save the output.
        append_to_file: If not None, then the data from the file named path_r
                     will also be appended to append_to_file
        var_name: Variable name used to save the data. Note that \"data\" is not a legal
                  name in stan.
        data_type: None to go by whether the experiment is binned, \"int\", or \"float\"."""
        if(data_type=="int" or
           (data_type is None and self.binned_expt)):
            tempdict = {var_name:self.data.astype(int)}
        else:
            tempdict = {var_name:self.data.astype(float)}
        try:
            stan_rdump(tempdict, path_r)
        except NameError as e:
            print("Error: %s; data could not be stored in an r file" % e)
        self._data_saved_loc = path_r
        if(not append_to_file is None):
            # Clean up this code
            file1 = open(path_r, 'r')
            file2 = open(append_to_file, 'r')
            text = file2.read() + "\n" + file1.read()
            file1.close()
            file2.close()
            outfile = open(append_to_file, 'w')
            outfile.write(text)
            outfile.clos()
        return

    def save_data_pdf(self, path_directory, path_pdf_prefix="data",
                      title="Data", y_label="Counts"):
        """ Save plots of the current data to pdf files
        For multidimensional data, it saves one row for each dimension, with the
        other dimensions set to the 0th index
        "_dimension_name.pdf" is appended to path_pdf_prefix
        The x label will be the name of the relevant dimension"""
        if not os.path.exists(path_directory):
            os.makedirs(path_directory)
        if(self.data is None):
            print("Data not set, and cannot be saved")
            return
        for i,d in enumerate(self.dimensions):
            curr_data = self._get_1d_data(i)
            x_data = curr_data[0,:]
            y_data = curr_data[1,:]
            x_label = d.name
            temp_gr = ROOT.TGraph(len(x_data),np.array(x_data),np.array(y_data))
            path_pdf = path_directory + "/"  + path_pdf_prefix + "_"
            path_pdf += d.name + ".pdf"
            rt.root_make_plot([temp_gr],path_pdf,["AP"],title,x_label,y_label,marker_styles=[3])
        self._data_saved_loc = path_directory
        return


    def log_likelihood_gauss_factor(self, param_mags=None):
        if(param_mags is None):
            param_mags = np.zeros(len(self.parameters))
            for i,p in enumerate(self.parameters):
                param_mags[i] = p.magnitude
        log_like = 0.0
        for i,p in enumerate(self.parameters):
            log_like += p.log_likelihood_gauss_factor(param_mags[i])
        return log_like

    def log_likelihood_poisson_factor(self, param_mags=None):
        if(param_mags is None):
            param_mags = np.zeros(len(self.parameters))
            for i,p in enumerate(self.parameters):
                param_mags[i] = p.magnitude
        exp_cts = self.calculate_expected_counts(param_mags)
        pois_log_probs = sp.stats.poisson.logpmf(self.data, exp_cts)
        return np.sum(pois_log_probs)

    def log_likelihood(self, param_mags=None):
        """Use the current data to minimize the data using either H1 or H0
        use_H0: Likelihood minimized with H0 if true, H1 otherwise
        param_mags = numpy array of the values of the parameters for which to calculate the magnitude
        Pre: data can be created via generate_fake_data(), or stored in data variable
        Post: """
        if(param_mags is None):
            param_mags = np.zeros(len(self.parameters))
            for i,p in enumerate(self.parameters):
                param_mags[i] = p.magnitude
        temp_mags = np.array(param_mags)
        gauss_log_like = self.log_likelihood_gauss_factor(param_mags)
        poisson_log_like = self.log_likelihood_poisson_factor(param_mags)      
        return gauss_log_like + poisson_log_like


class FrequentistAnalysis:
    """Class to perform a frequentist analysis"""

    def __init__(self, name, binned_analysis=True, experiments=list(),
                 parameters=list()):
        """Construct a new empty frequentist sensitivity analysis
        name: Name of the analysis, used when printing error messages
        binned_analysis: True if analysis should be binned, False if unbinned
                         Shapes and data must be regenerated if binned_analysis is modified
        The analysis is then initialized with the following methods:
          1. Use add_dimension to add all dimensions, such as time or energy
          2. Use add_parameter to add all parameters
          CONTINUE TO UPDATE THIS
        """
        self.name = name
        self.binned_analysis = binned_analysis
        self.experiments = experiments
        self.parameters = parameters
        return

    def __str__(self):
        str = "-----\n"
        str += "Analysis %s\n" % self.name
        str += "  All Parameters:\n"
        for p in self.parameters:
            str += "    %s\n" % p
        str += "  Experiments:\n"
        for e in self.experiments:
            str += "%s\n" % e
        str += "\n-----"
        return str

    def check_sanity(self):
        """Check that the parameter list for each experiment is the same
        as the parameter list for the FrequentistAnalysis object. Otherwise
        minimizing the log likelihood will not work as expected"""
        for e in self.experiments:
            if(not e.parameters==self.parameters):
                return False
        return True

    def get_parameter(self, name):
        """Get a Parameter object by name, so it can then be edited"""
        for p in self.parameters:
            if p.name==name:
                return p
                break
        else:
            print("Parameter \"%s\" not found by FrequentistAnalysis %s get_parameter()"%(name,self.name))  
        return

    def remove_parameter(self, name):
        """Remove the first Parameter object with the given name from all experiments
        and the FrequentistAnalyis object, and return it"""
        for i,p in enumerate(self.parameters):
            if p.name==name:
                for e in self.experiments:
                    e.remove_parameter(name)
                return self.parameters.pop(i)
                break
        else:
            print("Parameter \"%s\" not found by FrequentistAnalysis %s remove_parameter()"%(name,self.name))
            return None

    def print_parameters(self):
        """Print information about all parameters, including names used for access with get_parameter()"""
        str = "Parameters:\n"
        for p in self.parameters:
            str += "  %s\n" % p
        print(str)
        return

    
    def generate_all_shapes(self):
        """Generate shapes for parameters for all experiments"""
        for e in self.experiments:
            e.generate_all_shapes()
        return

    def save_all_shapes_root(self, path_root_prefix):
        """Save all shapes for all parameters in a rootfile with the given filename
        The tree for each shape will be the Shape name, the y branch will use the
        Parameter name, and the x branch will use the corresponding Dimension name
        (although any spaces replaced by "_" in the tree and branch names).
        Each Experiment will have its shapes saved in its own file, with file name
        path_root_prefix + "_" + experiment.name + ".root".
        """
        for e in self.experiments:
            path_root = path_root_prefix + "_"
            path_root += e.name.replace(" ","_") + ".root"
            e.save_all_shapes_root(path_root)
        return

    def save_all_shapes_pdf(self, path_directory_prefix):
        """Save pdf plots of all shapes.
        The filename and title will be the Shape name. The y label will use the
        Parameter name, and the x label will use the corresponding Dimension name.
        The shapes for each experiment will be placed in their own directory, with
        the directory name path_directory_prefix + "_" + experiment.name"""
        for e in self.experiments:
            path_directory = path_directory_prefix + "_" + e.name.replace(" ","_")
            e.save_all_shapes_pdf(path_directory)
        return
    
    def generate_fake_data(self):
        """Generate fake data for all experiments"""
        for e in self.experiments:
            e.generate_fake_data()
        return

    def save_all_data_root(self, path_root, br_name="counts"):
        """Save all current data in a root file
        The data for each experiment will be saved in its own tree.
        The tree name will be the experiment name."""
        for (i,e) in enumerate(self.experiments):
            if(i==0):
                e.save_data_root(path_root, tree_name_prefix=e.name.replace(" ","_"),
                                 br_name=br_name, recreate=True)
            else:
                e.save_data_root(path_root, tree_name_prefix=e.name.replace(" ","_"),
                                 br_name=br_name, recreate=False)
        return

    def save_all_data_r(self, path_r_prefix, var_name_prefix="counts",
                           append_to_file=None):
        """Save all current data in the r file format.
        The variable name for each experiment will be
        var_name_prefix + "_" + experiment.name
        The path name will be path_r_prefix + "_" + experiment.name + ".out\"
        If append_to_file is given, then all data will also be appended to the
        file named append_to_file"""
        for e in self.experiments:
            path_r = path_r_prefix + "_" + e.name.replace(" ","_") + ".out"
            var_name = var_name_prefix + "_" + e.name.replace(" ","_")
            e.save_data_r(path_r, var_name, append_to_file)
        return

    def save_all_data_pdf(self, path_directory, path_pdf_prefix="data",
                             title_prefix="Data", y_label="Counts"):
        """Save all current data in the given directory
        The experiment name will be attached to prefixes to make paths and titles"""
        for e in self.experiments:
            path_pdf_prefix2 = path_pdf_prefix + "_" + e.name.replace(" ","_")
            title = title_prefix + " " + e.name
            e.save_data_pdf(path_directory, path_pdf_prefix2, title, y_label)
        return


    def log_likelihood_total(self, param_mags=None):
        """Total log likehood for all experiments.
        if self.parameters has length N, then param_mags must have length N.
        generate_parameter_info() must be called first.
        If param_mags=None, then the magnitude of each parameter object is used
        to calculate the likelihood."""
        if(param_mags is None):
            param_mags = np.zeros(len(self.parameters))
            for i,p in enumerate(self.parameters):
                param_mags[i] = p.magnitude
        # Currently nothing is done to check that all parameters are inside their bounds
        # When I set the log_likelihood to 0 outside the bound, the minimizer had trouble
        # because outside the bounds there was nothing to tell it where to go.
        '''if(not len(param_mags)==len(self.parameters) or
           np.any(param_mags<self._param_lbs) or
           np.any(param_mags>self._param_ubs)):
            return 0'''
        log_like_tot = 0.0
        for e in self.experiments:
            log_like_tot += e.log_likelihood(param_mags)
        return log_like_tot

    def max_log_likelihood_fixed_params(self, param_indices_fixed, param_mags_fixed):
        """Total log likelihood of all experiments, in the case that some parameters are fixed
        and the likelihood is minimized over all other parameters.
        param_indices_fixed and param_mags_fixed the same length
        Post: returns a tuple containing an array with the values of all parameters (both those that 
                 were fixed and those that were minimized) in the first position, and the value of 
                 the max likelihood in the second position.
              A MinimizationFailedError may be thrown.
        """
        param_mags = np.zeros(len(self.parameters))
        for i,j in enumerate(param_indices_fixed):
            param_mags[j] = param_mags_fixed[i] # These will now remain fixed for all minimization
        param_indices_free = np.arange(0,len(self.parameters))
        param_indices_free = np.delete(param_indices_free, param_indices_fixed)
        def minus_ll_fixed(param_mags_free):
            for i,j in enumerate(param_indices_free):
                param_mags[j] = param_mags_free[i]
            return -self.log_likelihood_total(param_mags)
        x0 = np.zeros(len(param_indices_free))
        for i,j in enumerate(param_indices_free):
            x0[i] = self.parameters[j].magnitude
        res = sp.optimize.minimize(minus_ll_fixed, x0, method='nelder-mead')
        if(res.success):
            return (param_mags, -res.fun)
        else:
            raise MinimizaationFailedError("Minimization unsuccessful, output: %s" % res)

    def maximize_log_likelihood(self, use_H0=False, profile_nuisance=False):
        """Use the current data to minimize the data using either H1 or H0
        use_H0: Likelihood minimized with some parameters fixed at their H0 values
        profile_nuisance: Whether nuisance parameters should be profiled. If use_H0==False,
                          then all parameters are free and the global maximum likelihood
                          will be returned whether or not profile_nuisance is True
        Pre: data can be created via generate_fake_data(), or stored in data variable
        Post: returns a tuple containing an array with the values of all parameters that 
                give the maximum in the first position, and the value of 
                the max log_likelihood in the second position.
              A MinimizationFailedError may be thrown"""
        if(use_H0):
            if(profile_nuisance):
                H0_indices = list()
                H0_mags = list()
                non_nu_indices = list()
                x0 = list()
                for i,p in enumerate(self.parameters):
                    if(not p.magnitude_H0 is None):
                        H0_indices.append(i)
                        H0_mags.append(p.magnitude_H0)
                    elif(not p.nuisance):
                        non_nu_indices.append(i)
                        x0.append(p.magnitude)
                H0_indices = np.array(H0_indices)
                H0_mags = np.array(H0_mags)
                non_nu_indices = np.array(non_nu_indices)
                if(len(H0_indices)>0):
                    if(len(non_nu_indices)>0):
                        H0_non_nu_indices = np.append(H0_indices, non_nu_indices)
                    else:
                        H0_non_nu_indices = np.array(H0_indices)
                else:
                    H0_non_nu_indices = np.array(non_nu_indices)
                x0 = np.array(x0)
                
                all_params = [None]
                def minus_ll_prof(non_H0_non_nu_mags):
                    # minimize -log_likelihood with H0 and non_nuisance parameters fixed
                    if(len(H0_mags)>0):
                        if(len(non_H0_non_nu_mags)>0):
                            H0_non_nu_mags = np.append(H0_mags, non_H0_non_nu_mags)
                        else:
                            H0_non_nu_mags = np.array(H0_mags)
                    else:
                        H0_non_nu_mags = np.array(non_H0_non_nu_mags)
                    res = self.max_log_likelihood_fixed_params(H0_non_nu_indices, H0_non_nu_mags)
                    all_params[0] = res[0]
                    return -res[1]
                print(all_params[0])
                
                # maximize log_likelihood with only H0 fixed, and non_nuisace params profiled
                sp_res = sp.optimize.minimize(minus_ll_prof, x0, method='nelder-mead')
                res = (all_params[0], -sp_res.fun)
            else:
                H0_indices = list()
                H0_mags = list()
                for i,p in enumerate(self.parameters):
                    if(not p.magnitude_H0 is None):
                        H0_indices.append(i)
                        H0_mags.append(p.magnitude_H0)
                H0_indices = np.array(H0_indices)
                H0_mags = np.array(H0_mags)
                res = self.max_log_likelihood_fixed_params(H0_indices, H0_mags)
        else:
            # If H0 is not limiting parameters, then profiling does not matter
            res = self.max_log_likelihood_fixed_params(np.array([]),np.array([]))
        return res

    # Edit: For binned analysis, calculate chi2 for parameters to test goodness of fit
    def frequentist_analysis(self, path_root, iterations=100,
                             max_min_failures=100):
        my_file = ROOT.TFile(path_root, "RECREATE")
        res_tree = ROOT.TTree("results", "results")
        tmp_z = array('f', [0])
        res_tree.Branch("Z", tmp_z, "Z"+'/F')
        tmp_q0 = array('f', [0])
        res_tree.Branch("q0", tmp_q0, "q0"+'/F')
        params_tree = ROOT.TTree("params", "params")
        tmp_params = list()
        for i,p in enumerate(self.parameters):
            tmp_params.append(array('f', [0]))
            params_tree.Branch(p.name, tmp_params[i], p.name+'/F')

        self.generate_all_shapes()
        for it in range(0,iterations):
            self.generate_fake_data()
            # Profile nuisance should be true, but for very simple cases it doesn't matter, so
            # I can test like this
            ll_H0 = self.maximize_log_likelihood(use_H0=True, profile_nuisance=True)
            ll_H1 = self.maximize_log_likelihood(use_H0=False, profile_nuisance=False)
            # Edit: Use scipy.stats.chi2 pdf insted of assuing  chi2(1) distribution
            q0 = -2*(ll_H0[1]-ll_H1[1])
            tmp_q0[0] = q0
            #print(q0)
            Z = np.sqrt(q0)
            tmp_z[0] = Z
            res_tree.Fill()
            for i,arr in enumerate(tmp_params):
                tmp_params[i][0] = ll_H1[0][i]
            params_tree.Fill()
        res_tree.Write()
        params_tree.Write()
        my_file.Close()
        return
