# -*- coding: utf-8 -*-

"""extinction.py  - Classes for implementing detailed dust extinction models.

Old code, use with caution!
"""

import numpy as np
import warnings, sys
import scipy.interpolate as interp
from numpy.random import normal

thismod = sys.modules[__name__]


class Attenuator(object):
    """
    A level of abstraction around the dust curves, which allows for
    some generic methods. The Attenuator is composed of an (effective)
    DustLaw and a DustDistribution class, both of which can be easily
    extended.  The first takes wavelength and tau_v as arguments and
    returns \tau_\lambda. The second takes arguments that can be
    anything (e.g. stellar age, metallicity) and returns a
    distribution of tau_V
    """

    def __init__(self, dust_type='FM07'):
        self.name = dust_type
        self.dustcurve = getattr(thismod, dust_type)()
        # self.dustdist =

    def attenuate_spectrum(self, wave, inspec, pars):
        """This will be slow until the dust curves are vectorized.
        """
        spec = np.atleast_2d(inspec)  # memory hog
        for i, p in enumerate(pars):
            tau = self.dustcurve.acurve(wave, A_v=p['A_V'],
                                        R_v=p['R_V'],
                                        f_bump=p['F_BUMP'],
                                        uv_slope=p['UV_SLOPE'])
            spec[i, :] = spec[i, :]*np.exp(-tau)
        return spec

    def draw_taus(self, ntau, pars):
        """Not implemented
        """
        pass

# ------------------------------
# Extinction Law Classes
# ------------------------------


class GenericCurve(object):
    """Class for handling extinction curves with the parameterization
    of Fitzpatrick & Massa 1990 or F&M 2007
    """

    def __init__(self, form='F07', var=False):
        self.default_pars(var=var)

        if form == 'FM90':
            self._f = self.f_90
        elif form == 'F07':
            self._f = self.f_07

    def fm_curve(self, x, c1=None, c2=None, c3=None, c4=None,
                 c5=None, x0=4.59, gamma=0.89, **extras):
        """Fitzpatrick & Massa parameterization of the UV excess
        curve.
        """
        d = self.drude(x, gamma, x0)
        f = self._f(x, c5=c5)
        return (c1 + c2 * x + c3 * d + c4 * f)

    def drude(self, x, gamma, x0):
        """Drude profile for the 2175AA bump.

        :param x:
           Inverse wavelength (inverse microns) at which
           values for the drude profile are requested.

        :param gamma:
           Width of the Drude profile (inverse microns).

        :param x0:
           Center of the Drude profile (inverse microns).

        :returns k_lambda:
           The value of the Drude profile at x
        """
        return x**2 / ((x**2 - x0**2)**2 + (x * gamma)**2)

    def f_90(self, x, **extras):
        """UV rise in the excess curve, from FM90 (defined cubic
        polynomial).

        :param x:
            Inverse wavelength in inverse microns

        :returns:
            Cubic polynomial at x centered on 5.9 microns**(-1), as in
            FM90
        """
        return (x > 5.9)*(0.5392*(x-5.9)*(x-5.9)+0.05644*(x-5.9)**3)

    def f_07(self, x, c5=5.9, **extras):
        """UV rise in the excess curve, from FM07 (quadratic with
        variable center).

        :param x:
            Inverse wavelength in inverse microns.

        :param c5:
            Pivot point for the quadratic.

        :returns:
            Quadratic at x centered on c5, as in FM07.
        """
        return (x > c5)*(x-c5)*(x-c5)

    def powerlaw(self, x, R_v=3.1, k=None, alpha=-1.84, **extras):
        """Power-law shape of the color excess curve, for the NIR
        extinction curve, as in FM09.

        :param x:
            Inverse wavelength in inverse microns.

        :param R_v: (default: 3.1)
            Slope of the optical portion of the curve.

        :param k: (default: None)
            Power-law normalization.

        :param alpha: (default: -1.84)
            Power-law slope.
        """
        return k*(1./x)**alpha - R_v

    def spline(self, x, spline_x=[1.8, 2.5, 3.0],
               spline_k=[0.0, 1.32, 2.02]):
        """Cubic spline, for the optical portion of the extinction
        curve.

        :param x:
            Inverse wavelength in inverse microns.

        :param spline_x:
            Location of the spline anchor points, in inverse microns.

        :param spline_k:
            Value at the spline anchor points.
        """
        spline_x = np.asarray(spline_x)
        spline_k = np.asarray(spline_k)
        sp = interp.InterpolatedUnivariateSpline(spline_x, spline_k)
        return sp(x)

    def acurve(self, wave, **pars):
        """Return the extinction curve given the excess curve.
        """
        x = 1e4/wave
        self.setpars(**pars)
        e = self.ecurve(x, self.pardict)
        return e * self.pardict['A_v']/self.pardict['R_v'] + 1

    def setpars(self, A_v=1, R_v=3.1, f_bump=1.0, uv_slope=1,
                var=False, **extras):
        self.default_pars(var=var)
        self.pardict['c3'] *= f_bump
        self.pardict['A_v'] = A_v
        self.pardict['R_v'] = R_v
        self.pardict['c2'] *= uv_slope
        self.selfupdate(var=var)


class FM07(GenericCurve):
    """Extinction curves from Fitzpatrick and Massa 2007
    """

    def __init__(self, var=False):
        self.default_pars(var=var)
        self._f = self.f_07

    def default_pars(self, var=False):
        """FM07 parameterization coefficients, based on MW stars.

        :param var:
            Set to True to include scatter in the coefficients
            based on the FM07 sample.
        """
        p = {}
        # UV
        p['c2'] = 0.81
        p['c3'] = 2.99
        p['c4'] = 0.32
        p['c5'] = 6.10
        p['x0'] = 4.59
        p['gamma'] = 0.90
        # Optical
        p['ospline_x'] = 1e4/np.array([3300., 4000., 5330])
        p['ospline_k'] = (np.array([2.05, 1.32, 0.0]))
        # NIR
        p['alpha'] = 0-1.84
        p['R_v'] = 3.1

        if var:
            draws = normal(0, 1, 10)
            p['c2'] += draws[0] * 0.25
            p['c3'] += draws[1] * 1.0
            p['c4'] += draws[2] * 0.16
            p['c5'] += draws[3] * 0.6
            p['x0'] += draws[4] * 0.025
            p['gamma'] += draws[5] * 0.1
            p['ospline_k'] += draws[6:9] * np.array([0.17, 0.026, 0.008])
            p['R_v'] += draws[9] * 0.6

        self.pardict = p
        self.selfupdate(var=var)

    def selfupdate(self, var=False):
        """Enforce empirical relationships between some of the
        exctinction curve parameters, optionally including variance in
        the relationships.  These are based on data in FM09 and Gordon
        et al 2002.
        """
        k = -0.83 + 0.63 * self.pardict['R_v'] + var * normal(0, 0.12)
        c1 = 2.02 - 3.007 * self.pardict['c2'] + var * normal(0, 0.3)
        self.pardict['k'] = k
        self.pardict['c1'] = c1

    def ecurve(self, x, pardict):
        """Excess curve, including 3 components (UV, Optical, and
        NIR).
        """
        nir = self.powerlaw(x, **pardict)
        uv = self.fm_curve(x, **pardict)
        spline_x = np.array([1e4/2600, 1e4/2700] +
                            pardict['ospline_x'].tolist() +
                            [1.0, 0.75])
        spline_k = np.array(self.fm_curve(spline_x[0:2], **pardict).tolist() +
                            pardict['ospline_k'].tolist() +
                            self.powerlaw(spline_x[-2:], **pardict).tolist())

        optical = self.spline(x, spline_x=spline_x[::-1],
                              spline_k=spline_k[::-1])
        return (nir * (x < 0.75) + uv * (x > 3.8) +
                optical * ((x >= 0.75) & (x <= 3.8)))


class LMC(FM07):

    def __init__(self, var=False):
        self.default_pars(var=var)
        self._f = self.f_90

    def default_pars(self, var=False):
        """Gordon 2003 parameterization coefficients.  Set var to True
        to include scatter in the coefficients based on the Gordon 03
        sample.
        """
        p = {}
        # UV
        p['c2'] = 0.998 + var * normal(0, 0.027)
        p['c3'] = 2.719 + var * normal(0, 0.137)
        p['c4'] = 0.40 + var * normal(0, 0.036)
        p['x0'] = 4.579 + var * normal(0, 0.007)
        p['gamma'] = 0.934 + var * normal(0, 0.016)
        # Optical
        p['ospline_x'] = 1e4 / np.array([3300., 4000., 5330])
        p['ospline_k'] = (np.array([2.05, 1.32, 0.0]) +
                          var * normal(0, 1, 3) *
                          np.array([0.17, 0.026, 0.008]))
        # NIR
        p['alpha'] = -1.84
        p['R_v'] = 3.1 + var * normal(0, 0.6)

        self.pardict = p
        self.selfupdate()


class SMC(FM07):

    def __init__(self, var=False):
        self.default_pars(var=var)
        self._f = self.f_90

    def default_pars(self, var=False):
        """Gordon 2003 parameterization coefficients.  Set var to True
        to include scatter in the coefficients based on the Gordon 03
        sample.
        """
        p = {}
        # UV
        p['c2'] = 2.264 + var * normal(0, 0.040)
        p['c3'] = 0.389 + var * normal(0, 0.110)
        p['c4'] = 0.46 + var * normal(0, 0.079)
        p['x0'] = 4.60
        p['gamma'] = 1.0
        # Optical
        p['ospline_x'] = 1e4 / np.array([3300., 4000., 5330])
        p['ospline_k'] = (np.array([2.05, 1.32, 0.0]) +
                          var * normal(0, 1, 3) *
                          np.array([0.17, 0.026, 0.008]))
        # NIR
        p['alpha'] = -1.84
        p['R_v'] = 3.1 + var * normal(0, 0.6)

        self.pardict = p
        self.selfupdate()


class F99(GenericCurve):
    """Fitzpatrick 1999 R_v dependent extinction curves.  These are a
    one parameter family, though in practice we allow the bump
    strength to vary.
    """

    def __init__(self, var=False):
        self.default_pars(var=var)
        self._f = self.f_90

    def default_pars(self, var=False):
        """FM99 parameterization coefficients, for an R_v dependent
        curve.  Set var to True to include scatter in the coefficients
        based on the FM07 sample.
        """
        p = {}
        # NIR
        p['R_v'] = 3.1 + var * normal(0, 0.6)
        # Optical
        p['ospline_x'] = 1e4 / np.array([4100., 4670., 5470., 6000.,
                                         12200, 26500, 3e8])
        p['ospline_x'][-1] = 0.

        # UV
        p['c3'] = 3.23 + var * normal(0, 1.0)
        p['c4'] = 0.41 + var * normal(0, 0.16)
        p['x0'] = 4.596 + var * normal(0, 0.025)
        p['gamma'] = 0.99 + var * normal(0, 0.1)

        self.pardict = p
        self.selfupdate()

    def selfupdate(self, var=False):
        """Enforce empirical relationships between some of the
        exctinction curve parameters.
        """
        p = {}
        p['R_v'] = self.pardict['R_v']
        self.pardict['c2'] = (-0.824 + 4.717 / p['R_v'] +
                              var*normal(0, 0.25))
        self.pardict['c1'] = (2.02 - 3.007 * self.pardict['c2'] +
                              var * normal(0, 0.3))
        sk = np.array([1.208 + 0.0032 * p['R_v'] - 0.00033 * p['R_v']**2,
                       0.701 + 0.0016 * p['R_v'],
                       -0.050 + 0.0016 * p['R_v'],
                       -0.426 + 0.0044 * p['R_v'],
                       (0.829 / 3.1 - 1) * p['R_v'],
                       (0.265 / 3.1 - 1) * p['R_v'],
                       -p['R_v']])
        self.pardict['ospline_k'] = sk

    def ecurve(self, x, pardict):
        """Excess curve, including 2 components (UV and Optical/NIR)
        """
        uv = self.fm_curve(x, **pardict) * (x > 3.8)
        spline_x = np.array([1e4 / 2600, 1e4 / 2700] +
                            pardict['ospline_x'].tolist())
        spline_k = np.array(self.fm_curve(spline_x[0:2], **pardict).tolist() +
                            pardict['ospline_k'].tolist())
        optical = self.spline(x, spline_x=spline_x[::-1],
                              spline_k=spline_k[::-1]) * (x <= 3.8)
        return uv + optical
