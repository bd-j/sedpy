import numpy as np
import warnings, sys
import scipy.interpolate as interp


#### Attenuator class

# A level of abstraction around the dust curves, which allows for some
# generic methods.
# The Attenuator is composed of an (effective) DustLaw and a DustDistribution
# class, both of which can be easily extended.  The first takes wavelength and
# tau_v as arguments and returns tau_\lambda.  The second takes arguments
# that can be anything (e.g. stellar age, metallicity) and returns a
# distribution of tau_V

thismod = sys.modules[__name__]

class Attenuator(object):

    def __init__(self, dust_type = 'MilkyWay'):
        self.name = dust_type
        self.dustcurve = getattr(thismod, dust_type)()
        #self.dustdist = 

    def attenuate_spectrum(self, wave, inspec, pars):
        """This will be slow until the dust curves are vectorized"""
        spec = np.atleast_2d(inspec) #memory hog
        for i,p in enumerate(pars):
            tau = self.dustcurve.acurve(wave, A_v = p['A_V'], R_v = p['R_V'], f_bump = p['F_BUMP'])
            #print(tau.shape, spec.shape)
            spec[i,:] = spec[i,:]*np.exp(-tau)
        return spec

    def draw_taus(self,ntau, pars):
        pass
 

#########Extinction Law Classes
########

class Generic_Curve(object):
    """Class for handling extinction curves with the parameterization
    of Fitzpatrick & Massa 1990 or F&M 2007"""


    def __init__(self, form = 'F90', var = False):
        self.default_pars(var = var)
        if form is 'FM90':
            self._F = self.F_90
        else:
            self._F = self.F_07                        
 
    def fm_curve(self, x, c1 = None, c2 = None, c3 = None, c4 = None, c5 = None, x0 = 4.59, gamma  = 0.89):
        return  (c1 + c2*x + c3*self.drude(x, gamma, x0) + c4 * self._F(x, c5 = c5))
    
    def drude(self, x, gamma, x0):
        """Drude profile for the 2175AA bump"""
        return x*x / ( (x*x-x0*x0)**2 +(x*gamma)**2 )

    def F_90(self, x, **extras):
        """UV rise in the extinction curve.  defined cubic polynomial"""
        return (x > 5.9)*(0.5392*(x-5.9)*(x-5.9)+0.05644*(x-5.9)**3)

    def F_07(self, x, c5 = 5.9, **extras):
        """UV rise in the extinction curve, FM07 (simple quadratic with variable center)"""
        return (x > c5)*(x-c5)*(x-c5)

    def powerlaw(self, x, R_v = 3.1, k = None, alpha = -1.84):
        return k*(1./x)**alpha - R_v

    def spline(self, x, spline_x = [1.8,2.5,3.0], spline_k = [0.0,1.32, 2.02]):
        spline_x = np.asarray(spline_x)
        spline_k = np.asarray(spline_k)
        sp = interp.InterpolatedUnivariateSpline(spline_x, spline_k)
        return sp(x)
        
    def acurve(self, wave, **pars):
        x = 1e4/wave
        self.setpars(pars)
        return self.ecurve(x, self.pardict) * self.pardict['A_v']/self.pardict['R_v'] + 1


class F99(Generic_Curve):
    pass

class LMC(Generic_Curve):
    pass

class SMC(Generic_Curve):
    pass

class FM07(Generic_Curve):
    """ Extinction curves from Fitzpatrick and Massa 2007"""

    def default_pars(self, var = False):
        """FM99 parameterization coefficients."""
        from numpy.random import normal 
        p = {}
        #UV
        p['c2'] = 0.81 + var*normal(0, 0.25)
        p['c1'] = 2.02 - 3.007 * p['c2'] +var*normal(0, 0.3)
        p['c3'] = 2.99 + var*normal(0, 1.0)
        p['c4'] = 0.32 +var*normal(0, 0.16)
        p['c5'] = 6.10 + var*normal(0, 0.6)
        p['x0'] = 4.59+ var*normal(0, 0.025)
        p['gamma'] = 0.90 + var*normal(0, 0.1)
        #Optical
        p['ospline_x'] = 1e4/np.array([3300., 4000., 5330])
        p['ospline_k'] = (np.array([2.05, 1.32, 0.0])
                          + var*normal(0,1,3)*np.array([0.17, 0.026, 0.008]))
        #NIR
        p['alpha'] = -1.84
        p['R_v'] = 3.1 +var*normal(0., 0.6)
        p['k'] = -0.83 + 0.63*p['R_v'] + var*normal(0, 0.12)
        
        self.pardict = p

    def selfupdate(self, var = False):
        self.pardict['k'] = -0.83 + 0.63*self.pardict['R_v'] + var*normal(0, 0.12)
        self.pardict['c1'] = 2.02 - 3.007 * self.pardict['c2'] +var*normal(0, 0.3)
        
    def setpars(self, inpars):
        self.pardict['c3'] *= inpars['f_bump']
        self.pardict['A_v'] = inpars['A_v']
        self.pardict['R_v'] = inpars['R_v']
        
        self.selfupdate()
        
    def ecurve(self, x, pardict):
        nir = self.powerlaw( x, **pardict) * (x < 0.75)
        uv = self.fm_curve(x, **pardict) * (x > 3.8)
        spline_x = np.array([1e4/2600, 1e4/2700, pardict['ospline_x'], 1.0,0.75])
        spline_k = [self.fm_curve(spline_x[0:2], **pardict),
                    pardict['ospline_k'], self.powerlaw( spline_x[-2:], **pardict)]
        optical  = self.spline(x, spline_x =spline_x, spline_k = spline_k) * ( (x >= 0.75) & (x <= 3.8)) 
        return nir + uv + optical

    
        


