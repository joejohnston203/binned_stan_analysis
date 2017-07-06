# CNSexperiment.py
#
# Class that represents a CNS experiment with some tools for generating fake
# data and computing limits.
#
# Adam Anderson
# 12 April 2016
# adama@fnal.gov
#
# Note: Convention on units:
#   --all masses are in kg
#   --all energies are in keV
#   --all distances are in cm

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
import scipy.integrate as spint
import pdb
import ReactorTools

# fundamental constants
keVPerGeV       = 1e6      # [keV / GeV]
keVPerMeV       = 1e3
gPerkg          = 1e3
hbarc 		= 0.197*keVPerGeV	# [keV fm]
fmPercm		= 1.0e13	# [fm / cm]
sin2thetaW	= 0.2387
cLight		= 3.0e8 	# [m / s]
nAvogadro	= 6.022e23
joulePereV	= 1.602e-19	# [J / eV]
eVPerFission	= 200.0e6 	# [eV]
Mn              = 0.931 * keVPerGeV
Gfermi		= (1.16637e-5 / (keVPerGeV**2.))*(hbarc/fmPercm)	# [cm / keV]


class CNSexperiment:
    def __init__(self, N, Z, dRdEnu, T_bg, dRdT_bg, detMass, time, nuFlux):
        '''
        Constructor for CNSexperiment object, which handles rate calculations,
        toy MC, and fitting.

        Parameters
        ----------
        N : float
            Number of neutrons
        Z : float
            Number of protons
        dRdEnu : func
            Differential neutrino spectrum per fission, as a function of
            neutrino energy
        T_bg : array
            Recoil energies corresponding to rates in dRdT_bg; background rates
            are interpolated from sampled points, and this is energy of each
            sample
        dRdT_bg : function
            Differential background rate per recoil energy, as a function of
            recoil energy; background rates are interpolated from these sampled
            points
        detMass : float
            Detector mass in kg
        time : float
            Live time in seconds
        nuFlux : float
            Total neutrino flux at the detector in nu / cm^2 / s

        Returns
        -------
        None
        '''
        self.nNeutrons 		= N
        self.nProtons 		= Z
        self.nNucleons          = N + Z
        self.isotopeMassc2	= self.nNucleons * Mn 		# [keV]
        self.dRdEnu_source      = dRdEnu
        self.T_background       = T_bg
        self.dRdT_background    = dRdT_bg
        self.Qweak		= self.nNeutrons - (1.0 - 4.0*sin2thetaW) * self.nProtons
        self.detectorMass       = detMass
        self.nTargetAtoms       = self.detectorMass * gPerkg / self.nNucleons * nAvogadro # approximately
        self.livetime           = time
        self.nuFlux             = nuFlux

        # evaluate the cdf at test points and invert by linear interpolation
        self.T_for_cdf = np.insert(np.logspace(-3,3,100), 0, 1e-8)
        self.dRdT_samples = self.dRdT_CNS(self.T_for_cdf)
        self.cumul_rate = np.array([spint.trapz(self.dRdT_samples[:jT], self.T_for_cdf[:jT]) for jT in np.arange(1,len(self.T_for_cdf))])
        self.total_rate = self.cumul_rate[-1]
        self.cdf = self.cumul_rate / self.total_rate


    def dsigmadT_atEnu_CNS(self, Enu, T):
        '''
        Differential cross section per recoil energy as a function of recoil
        energy at a **fixed** neutrino energy (monoenergetic neutrinos) for CNS.

        Parameters
        ----------
        Enu : float or array
            Neutrino energy (keV)
        T : float
            Recoil energy (keV)

        Returns
        -------
        dsigmadT : array
            Differential cross section (cm^2/keV)
        '''
        # ensure that Enu is always a numpy array, even if input type is float or python list
        if type(Enu) == float:
    		Enu = np.asarray([Enu])
    	else:
    		Enu = np.asarray(Enu)

        dsigmadT = Gfermi**2. / (4*np.pi) * self.Qweak**2. * self.isotopeMassc2 * \
                   (1. - (self.isotopeMassc2 * T) / (2.*Enu**2.)) * self.F_Helm(T, self.nNucleons)
#        dsigmadT = (2.30156e-25)**2. / (4*np.pi) * (39.1836)**2. * (6.765e7) * \
#                    (1. - ((6.765e7)* T) / (2.*Enu**2.)) * 1.
        dsigmadT[dsigmadT<0] = 0.
#        print(T)
#        print(Enu)
#        print(Gfermi)
#        print(self.Qweak)
#        print(self.isotopeMassc2)
#        print((1. - (self.isotopeMassc2 * T) / (2.*Enu**2.)))
#        print(self.F_Helm(T, self.nNucleons))
        return dsigmadT


    def dsigmadT_CNS(self, T):
        '''
        Differential cross section, integrated over the energy spectrum of
        incident neutrinos. Note that scipy's numerical integrator "quad" is
        quite slow, so this function does not evaluate quickly (not good for
        use in optimization algorithms).

        Parameters
        ----------
        T : float or array
            Recoil energy (keV)

        Returns
        -------
        dsigmadT : array
            Differential cross section 
        '''
        # ensure that T is always a numpy array, even if input type is float or python list
        if type(T) == float:
    		T = np.asarray([T])
    	else:
    		T = np.asarray(T)
        def dsigmadTdEnu_CNS(Enu, T):
            dsigmadTdEnu = 1.e42 * self.dsigmadT_atEnu_CNS(Enu, T) * self.dRdEnu_source(Enu)
            return dsigmadTdEnu

        # pdb.set_trace()
        dsigmadT = np.array([spint.quad(dsigmadTdEnu_CNS, 0, 1.e6, args=(Tval), epsabs=1.e-6)[0]/1.e42 for Tval in T])
        return dsigmadT


    def dRdT_CNS(self, T):
        '''
        Differential rate per recoil energy. This calls functions that use scipy
        "quad", which is slow (see alternative "dRdT_CNS_fast" below).

        Parameters
        ----------
        T : float or array
            Recoil energy (keV)

        Returns
        -------
        dRdT : array
            Differential rate (evts/keV/kg)
            (we multiply by time, so this is the differential rate for the entire runtime)
        '''
        dRdT = self.dsigmadT_CNS(T) * self.nuFlux * self.livetime * self.nTargetAtoms
        return dRdT


    def dRdT_CNS_fast(self, T):
        '''
        Differential rate per recoil energy from interpolating function. We
        evaluate dRdT by interpolating the more accurate version above.
        Interpolation is only valid between 1e-8 keV and 100 keV.

        Parameters
        ----------
        T : float or array
            Recoil energy

        Returns
        -------
        dRdT : array
            Differential rate
        '''
        dRdT = np.interp(T, self.T_for_cdf, self.dRdT_samples)
        return dRdT


    def dRdT_bg(self, T):
        '''
        Differential background rate.

        Parameters
        ----------
        T : float or array
            Recoil energy

        Returns
        -------
        dRdT : array
            Differential rate
        '''
        dRdT = np.interp(T, self.T_background, self.dRdT_background)
        return dRdT


    def F_Helm(self, T, A):
        '''
        Helm form factor (same as Lewin and Smith DM paper).

        Parameters
        ----------
        T : float or array
            Recoil energy
        A : float
            Number of nucleons

        Returns
        -------
        F : float or array
            Form factor
        '''
        # define the momentum transfer in MeV / c
        q = np.sqrt(2 * A * Mn * T)

        # Standard Helm form factor
        R = 0.89*A**(1.0/3.0) + 0.3   # nuclear radius [fm]
        a = 0.52                      # [fm]
        s = 0.9                       # smearing parameter [fm]
        c = 1.23 * A**(1.0/3.0) - 0.6 # [fm]
        rn = np.sqrt(c**2.0 + (7.0/3.0)*(np.pi*a)**2 - 5*s**2)   # [fm]
        F = (3 / ((q * rn / hbarc)**3) * (np.sin(q * rn / hbarc) - (q * rn / hbarc) * np.cos(q * rn / hbarc)))**2.0 * \
            np.exp(-1.0 * (q*s / hbarc)**2)

        return F


    def neg2logL(self, mu, MCdata):
        '''
        The -2logL as a function of a signal strength parameter mu. The
        parameter mu is defined such that mu=1 corresponds to the nominal cross
        section of the CNS signal.

        Parameters
        ----------
        mu : float
            Signal strength parameter

        Returns
        -------
        n2LogL : float
        '''
        bg_rate = 0.
        sig_rate = mu*self.total_rate

        n2LogL = 2.*(bg_rate + sig_rate) - \
                 2.*np.sum(np.log(sig_rate * self.dRdT_CNS_fast(MCdata) + \
                                  bg_rate * self.dRdT_bg(MCdata)))

        return n2LogL


    def run_toy(self):
        '''
        Generates a pseudoexperiment by toy MC.

        Parameters
        ----------
        None

        Returns
        -------
        T : array
            Recoil energy of MC samples
        '''
        # use inverse CDF method for random sampling of the energy spectrum
        nEvt = np.random.poisson(lam=self.total_rate)
        unirand = np.random.random(nEvt)

        # evaluate inverse of cdf by linear interpolation
        T = np.interp(unirand, self.cdf, self.T_for_cdf[:-1])

        return T
