# ReactorTools.py
#
# Some tools for calculating the neutrino rate from a nuclear reactor.
#
# Adam Anderson
# 14 April 2016
# adama@fnal.gov
#
# Note: Convention on units:
#   --all masses are in kg
#   --all energies are in keV

import numpy as np
import ROOT

def dRdEnu_U235(Enu):
	'''
	Reactor anti neutrino spectrum from U235 (see arXiv:1101.2663v3), per
	fission

	Parameters
	----------
	Enu : array
		Neutrino energy in keV

	Returns
	-------
	spectrum : array
		Spectrum [nu / keV / fission]
	'''
	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
	EnuMeV = Enu / 1.e3
	spectrum = 1e-3 * np.exp(3.217 - 3.111*EnuMeV + 1.395*(EnuMeV**2.0) - \
					  (3.690e-1)*(EnuMeV**3.0) + (4.445e-2)*(EnuMeV**4.0) - (2.053e-3)*(EnuMeV**5.0))
	#spectrum[EnuMeV<1.0] = 1e-3 * np.exp(3.217 - 3.111*1.0 + 1.395*(1.0**2.0) - \
#					  (3.690e-1)*(1.0**3.0) + (4.445e-2)*(1.0**4.0) - (2.053e-3)*(1.0**5.0))
	return spectrum


def dRdEnu_U238(Enu):
	'''
	Reactor anti neutrino spectrum from U238 (see arXiv:1101.2663v3), per
	fission

	Parameters
	----------
	Enu : array
		Neutrino energy in keV

	Returns
	-------
	spectrum : array
		Spectrum [nu / keV / fission]
	'''
	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
	EnuMeV = Enu / 1.e3
	spectrum = 1e-3 * np.exp((4.833e-1) + (1.927e-1)*EnuMeV - (1.283e-1)*EnuMeV**2.0 - \
						(6.762e-3)*EnuMeV**3.0 + (2.233e-3)*EnuMeV**4.0 - (1.536e-4)*EnuMeV**5.0)
	#spectrum[EnuMeV<1.0] = 1e-3 * np.exp((4.833e-1) + (1.927e-1)*1.0 - (1.283e-1)*1.0**2.0 - \
	#					(6.762e-3)*1.0**3.0 + (2.233e-3)*1.0**4.0 - (1.536e-4)*1.0**5.0)
	return spectrum


def dRdEnu_Pu239(Enu):
	'''
	Reactor anti neutrino spectrum from Pu239 (see arXiv:1101.2663v3), per
	fission

	Parameters
	----------
	Enu : array
		Neutrino energy in keV

	Returns
	-------
	spectrum : array
		Spectrum [nu / keV / fission]
	'''
	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
	EnuMeV = Enu / 1.e3
	spectrum = 1e-3 * np.exp(6.413 - 7.432*EnuMeV + 3.535*EnuMeV**2.0 - \
						(8.82e-1)*EnuMeV**3.0 + (1.025e-1)*EnuMeV**4.0 - (4.550e-3)*EnuMeV**5.0)
	#spectrum[EnuMeV<1.0] = 1e-3 * np.exp(6.413 - 7.432*1.0 + 3.535*1.0**2.0 - \
	#					(8.82e-1)*1.0**3.0 + (1.025e-1)*1.0**4.0 - (4.550e-3)*1.0**5.0)
	return spectrum

def dRdEnu_Pu241(Enu):
	'''
	Reactor anti neutrino spectrum from Pu239 (see arXiv:1101.2663v3), per
	fission

	Parameters
	----------
	Enu : array
		Neutrino energy in keV

	Returns
	-------
	spectrum : array
		Spectrum [nu / keV / fission]
	'''
	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
	EnuMeV = Enu / 1.e3
	spectrum = 1e-3 * np.exp(3.251 - 3.204*EnuMeV + 1.428*EnuMeV**2.0 - \
						(3.675e-1)*EnuMeV**3.0 + (4.254e-2)*EnuMeV**4.0 - (1.896e-3)*EnuMeV**5.0)
	#spectrum[EnuMeV<1.0] = 1e-3 * np.exp(3.251 - 3.204*1.0 + 1.428*1.0**2.0 - \
	#					(3.675e-1)*1.0**3.0 + (4.254e-2)*1.0**4.0 - (1.896e-3)*1.0**5.0)
	return spectrum

def dRdEnu_fuel_frac(Enu):
	'''
	Reactor anti neutrino spectrum from Pu239 (see arXiv:1101.2663v3), per
	fission

	Parameters
	----------
	Enu : array
		Neutrino energy in keV

        The fuel fractions are hard coded for now, but eventually I will
        figure out how to include these as parameters

	Returns
	-------
	spectrum : array
		Spectrum [nu / keV / fission]
	'''
	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
        U235frac = 0.556;
        U238frac = 0.071;
        Pu239frac = 0.326;
        Pu241frac = 0.047;

	spectrum = U235frac*dRdEnu_U235(Enu) + U238frac*dRdEnu_U238(Enu) +\
                   Pu239frac*dRdEnu_Pu239(Enu) + Pu241frac*dRdEnu_Pu241(Enu)

	return spectrum


def nuFlux(power, distance):
	'''
	Computes the total flux per fission of reactor antineutrinos at a given
	distance from the core, assuming a point-like flux, and nominal neutrino production

	Parameters
	----------
	power : float
		Reactor power in MW
	distance : float
		Distance in cm from reactor core at which flux is to be calculated

	Returns
	-------
	flux : float
		The reactor neutrino flux in fissions/s/cm^2 
	'''
	flux = power/200.0/1.602176565e-19 / (4*np.pi * distance**2.)
	return flux

# Setup for the huber spectra
huber_setup_complete = False
spl_U235 = 0
spl_Pu239 = 0
spl_Pu241 = 0
def spl_U235_eval(): return
def spl_Pu239_eval(): return
def spl_Pu241_eval(): return

def Huber_setup(file_U235='/home/joe/morpho-develop/examples/Ricochet/data/U235-anti-neutrino-flux-250keV.dat',
                file_Pu239="/home/joe/morpho-develop/examples/Ricochet/data/Pu239-anti-neutrino-flux-250keV.dat",
                file_Pu241="/home/joe/morpho-develop/examples/Ricochet/data/Pu241-anti-neutrino-flux-250keV.dat"):
        global huber_setup_complete
        global spl_U235
        global spl_Pu239, spl_Pu241
        global spl_U235_eval, spl_Pu239_eval, spl_Pu241_eval
        # U235
        enU235, specU235 = np.loadtxt(file_U235,usecols=(0,1),unpack=True)
        spl_U235 = ROOT.TGraph(len(enU235),np.ascontiguousarray(enU235),np.ascontiguousarray(specU235))
        def spl_U235_temp_eval(Enu):
                return spl_U235.Eval(Enu)
        spl_U235_eval = np.vectorize(spl_U235_temp_eval)
        # Pu239
        enPu239, specPu239 = np.loadtxt(file_Pu239,usecols=(0,1),unpack=True)
        spl_Pu239 = ROOT.TGraph(len(enPu239),np.ascontiguousarray(enPu239),np.ascontiguousarray(specPu239))
        def spl_Pu239_temp_eval(Enu):
                return spl_Pu239.Eval(Enu)
        spl_Pu239_eval = np.vectorize(spl_Pu239_temp_eval)
        # Pu241
        enPu241, specPu241 = np.loadtxt(file_Pu241,usecols=(0,1),unpack=True)
        spl_Pu241 = ROOT.TGraph(len(enPu241),np.ascontiguousarray(enPu241),np.ascontiguousarray(specPu241))
        def spl_Pu241_temp_eval(Enu):
                return spl_Pu241.Eval(Enu)
        spl_Pu241_eval = np.vectorize(spl_Pu241_temp_eval)
        # Finished
        huber_setup_complete = True

def dRdEnu_U235_Huber(Enu):
        # check global setup
        global huber_setup_complete
        global spl_U235_eval
        if(not huber_setup_complete):
                Huber_setup()

	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
        # input is in keV, huber spline expects MeV
        # huber spline gives results in 1/Mev/fission, we want 1/keV/fission
        spec = spl_U235_eval(Enu*1e-3)*1e-3
        spec[Enu<2.e3] = spl_U235_eval(2.0)*1e-3
        spec[spec<0] = 0
        return spec

def dRdEnu_Pu239_Huber(Enu):
        # check global setup
        global huber_setup_complete
        global spl_Pu239_eval
        if(not huber_setup_complete):
                Huber_setup()

	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
        # input is in keV, huber spline expects MeV
        # huber spline gives results in 1/Mev/fission, we want 1/keV/fission
        spec = spl_Pu239_eval(Enu*1e-3)*1e-3
        spec[Enu<2.e3] = spl_Pu239_eval(2.0)*1e-3
        spec[spec<0] = 0
        return spec

def dRdEnu_Pu241_Huber(Enu):
        # check global setup
        global huber_setup_complete
        global spl_Pu241_eval
        if(not huber_setup_complete):
                Huber_setup()

	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
        # input is in keV, huber spline expects MeV
        # huber spline gives results in 1/Mev/fission, we want 1/keV/fission
        spec = spl_Pu241_eval(Enu*1e-3)*1e-3
        spec[Enu<2.e3] = spl_Pu241_eval(2.0)*1e-3
        spec[spec<0] = 0
        return spec

def dRdEnu_U238_Huber(Enu):
        # Huber does not have data for U238
        return dRdEnu_U238(Enu)

def dRdEnu_fuel_frac_Huber(Enu):
	'''
	Reactor anti neutrino spectrum from Pu239 (see arXiv:1101.2663v3), per
	fission

	Parameters
	----------
	Enu : array
		Neutrino energy in keV

        The fuel fractions are hard coded for now, but eventually I will
        figure out how to include these as parameters

	Returns
	-------
	spectrum : array
		Spectrum [nu / keV / fission]
	'''
	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)

        U235frac = 0.556;
        U238frac = 0.071;
        Pu239frac = 0.326;
        Pu241frac = 0.047;

	spectrum = U235frac*dRdEnu_U235_Huber(Enu) +\
                   U238frac*dRdEnu_U238_Huber(Enu) +\
                   Pu239frac*dRdEnu_Pu239_Huber(Enu) +\
                   Pu241frac*dRdEnu_Pu241_Huber(Enu)

	return spectrum
