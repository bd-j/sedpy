# -*- coding: utf-8 -*-

"""attenuation.py  - A collection of simple functions implementing various
proposed attenuation curves.
"""

import numpy as np
import warnings, sys

# --------------------
# ATTENUATION CURVES
# --------------------

__all__ = ["calzetti", "chevallard", "conroy", "noll",
           "powerlaw", "drude", "broken_powerlaw",
           "cardelli", "smc", "lmc"]


def powerlaw(wave, tau_v=1, alpha=1.0, **kwargs):
    r"""Simple power-law attenuation, normalized to 5500\AA.

    :param wave:
        The wavelengths at which optical depth estimates are desired.

    :param tau_v: (default: 1)
        The optical depth at 5500\AA, used to normalize the
        attenuation curve.

    :returns tau:
        The optical depth at each wavelength.
    """
    return tau_v * (wave / 5500)**(-alpha)


def calzetti(wave, tau_v=1, R_v=4.05, **kwargs):
    r"""Calzetti et al. 2000 starburst attenuation curve, with extrapolations to
    the FUV and NIR.

    :param wave:
        The wavelengths at which optical depth estimates are desired.

    :param tau_v: (default: 1)
        The optical depth at 5500\AA, used to normalize the attenuation curve.

    :param R_v: (default: 4.05)
        The ratio of total selective extinction, parameterizing the slope of
        the attenuation curve.  A_v = R_v * E(B-V)

    :returns tau:
        The optical depth at each wavelength.
    """
    # optical/NIR
    k1 = lambda x: 2.659 * (-1.857 + 1.040 * x)
    # UV
    k2 = lambda x: 2.659 * (-2.156 + 1.509 * x - 0.198 * x**2. + 0.011 * x**3.)

    # get slopes at edges and k(5500)
    uv = np.array([0.12, 0.13]) * 1e4
    kuv = k2(1e4 / uv) + R_v
    uv_slope = np.diff(kuv) / np.diff(uv)
    ir = np.array([2.19, 2.20]) * 1e4
    kir = k1(1e4 / ir) + R_v
    ir_slope = np.diff(kir) / np.diff(ir)
    k_v = k2(1e4 / 5500.) + R_v

    # define segments
    uinds = (wave >= 1200.) & (wave < 6300)  # uv
    oinds = (wave >= 6300.) & (wave <= 22000)  # optical
    xinds = (wave < 1200.)  # xuv
    iinds = (wave > 22000.)  # ir

    # do it
    x = 1e4 / wave
    ktot = oinds * (k1(x) + R_v)
    ktot += uinds * (k2(x) + R_v)
    ktot += xinds * (kuv[0] + (wave - uv[0]) * uv_slope)
    ktot += iinds * (kir[1] + (wave - ir[1]) * ir_slope)

    ktot[ktot < 0] = 0
    tau_lambda = tau_v * (ktot / k_v)
    return tau_lambda


def drude(x, x0=4.59, gamma=0.90, **extras):
    """Drude profile for the 2175AA bump.

    :param x:
        Inverse wavelength (inverse microns) at which values for the drude
        profile are requested.

    :param gamma:
        Width of the Drude profile (inverse microns).

    :param x0:
       Center of the Drude profile (inverse microns).

    :returns k_lambda:
       The value of the Drude profile at x, normalized such that the peak is 1.
     """
    #return (w * gamma)**2 / ((w**2 - w0**2)**2 + (w * gamma)**2)
    return (x*gamma)**2 / ((x**2 - x0**2)**2 + (x * gamma)**2)


def noll(wave, tau_v=1, delta=0.0, c_r=0.0, Ebump=0.0, **kwargs):
    r"""Noll 2009 attenuation curve. This is based on the Calzetti curve, with
    added variable bump (as Drude) and overall slope change.  Any extra
    keywords are passed to the Drude (e.g. x0, gamma, both in inverse microns).

    :param wave:
        The wavelengths at which optical depth estimates are desired.

    :param tau_v: (default: 1)
        The optical depth at 5500\AA, used to normalize the attenuation curve.

    :param Ebump: (default: 0.0)
        Stength of the 2175\AA bump.  Normalizes the Drude profile.

    :param delta: (default 0.)
        Slope of the power-law that modifies the Calzetti curve.

    :param c_r:
        Constant used to alter R_v=A_V/E(B-V) of the curve.  To maintain the
        Calzetti R_v=4.05, use c_r = -delta.  Note that even with c_r = -delta
        the calzetti curve will not be recovered unless delta=0

    :returns tau:
        The optical depth at each wavelength.
    """
    kcalz = calzetti(wave, tau_v=1.0, R_v=4.05) - 1
    k = kcalz + Ebump / 4.05 * drude(1e4 / wave, **kwargs)
    a = (k * (1 - 1.12 * c_r) + 1) * (wave / 5500.)**delta
    return a * tau_v


def chevallard(wave, tau_v=1, **kwargs):
    r""" \tau_v dependent attenuation curves matched to disk RT models,
    as in Chevallard et al. 2013.  No UV bump (or indeed tests in the
    UV at all).

    :param wave:
        The wavelengths at which optical depth estimates are desired.

    :param tau_v: (default: 1)
        The optical depth at 5500\AA, used to normalize the
        attenuation curve.

    :returns tau:
        The optical depth at each wavelength.
    """
    # missing a UV bump
    alpha_v = 2.8 / (1 + np.sqrt(tau_v))  # +/- 25%
    bb = 0.3 - 0.05 * tau_v  # +/- 10%
    alpha = alpha_v + bb * (wave * 1e-4 - 0.55)
    tau_lambda = tau_v * (wave / 5500.0)**(-alpha)
    return tau_lambda


def conroy(wave, tau_v=1, R_v=3.1, f_bump=0.6, **kwargs):
    r""" Conroy & Schiminovich 2010 dust attenuation curves including a
    decreased UV bump.

    :param wave:
        The wavelengths at which optical depth estimates are desired.

    :param tau_v: (default: 1)
        The optical depth at 5500\AA, used to normalize the
        attenuation curve.

    :param R_v: (default: 3.1)
        The ratio of total selective extinction, parameterizing the
        slope of the attenuation curve.  A_v = R_v * E(B-V)

    :param f_bump: (default: 0.6)
        The strength of the 2175\AA UV bump, as a fraction of the bump
        strength in Cardelli et al. extinction curve.

    :returns tau:
        The optical depth at each wavelength.
    """
    x = 1e4 / wave
    nx = x.shape[0]
    a = np.zeros_like(x)
    b = np.zeros_like(x)

    # IR 0.909 - 3.3 micron
    ir = (x >= 0.3) & (x < 1.1)
    a[ir] = 0.574 * x[ir]**1.61
    b[ir] = -0.527 * x[ir]**1.61

    # optical 0.303 - 0.909 micron
    opt = (x >= 1.1) & (x < 3.3)
    y = x[opt]-1.82
    a[opt] = (1 + 0.177 * y - 0.504 * y**2 - 0.0243 * y**3 +
              0.721 * y**4 + 0.0198 * y**5 - 0.7750 * y**6 +
              0.330 * y**7)
    b[opt] = (1.413 * y + 2.283 * y**2 + 1.072 * y**3 -
              5.384 * y**4 - 0.622 * y**5 + 5.303 * y**6 -
              2.090 * y**7)

    # NUV 0.17 to 0.303 micron
    nuv = (x >= 3.3) & (x < 5.9)
    tmp = (-0.0370 + 0.0469 * f_bump - 0.601 * f_bump / R_v + 0.542 / R_v)
    fa = (3.3 / x[nuv])**6. * tmp
    tmp = 0.104 * f_bump / ((x[nuv] - 4.67)**2 + 0.341)
    a[nuv] = 1.752 - 0.316 * x[nuv] - tmp + fa
    tmp = 1.206 * f_bump / ((x[nuv] - 4.62)**2 + 0.263)
    b[nuv] = -3.09 + 1.825 * x[nuv] + tmp

    # FUV 0.125 - 0.17 micron
    fuv = (x >= 5.9) & (x < 8.0)
    fa = -0.0447 * (x[fuv] - 5.9)**2.0 - 0.00978 * (x[fuv] - 5.9)**3
    fb = 0.213 * (x[fuv] - 5.9)**2. + 0.121 * (x[fuv] - 5.9)**3
    tmp = 0.104 * f_bump / ((x[fuv] - 4.67)**2 + 0.341)
    a[fuv] = 1.752 - 0.316 * x[fuv] - tmp + fa
    tmp = 1.206 * f_bump / ((x[fuv] - 4.62)**2 + 0.263)
    b[fuv] = -3.09 + 1.825 * x[fuv] + tmp + fb

    alam = (a + b / R_v)

    # XUV below 1250AA
    xuv = x >= 8.0
    x8 = 8.0
    fa = -0.0447 * (x8 - 5.9)**2 - 0.00978 * (x8 - 5.9)**3
    fb = 0.213 * (x8 - 5.9)**2. + 0.121 * (x8 - 5.9)**3
    tmp = 0.104 * f_bump / ((x8 - 4.67)**2 + 0.341)
    af = 1.752 - 0.316 * x8 - tmp + fa
    tmp = 1.206 * f_bump / ((x8 - 4.62)**2 + 0.263)
    bf = -3.09 + 1.825 * x8 + tmp + fb
    a8 = (af + bf / R_v)

    alam[xuv] = (x8 / x[xuv])**(-1.3) * a8

    return tau_v * alam


def broken_powerlaw(wave, tau_v=1, alpha=[0.7, 0.7, 0.7],
                    breaks=[0, 3000, 10000, 4e4], **kwargs):
    r""" Attenuation curve as in V. Wild et al. 2011, i.e. power-law
    slope can change between regions.  Superceded by Chevallard 2013
    for optical/NIR.

    :param wave:
        The wavelengths at which optical depth estimates are desired.

    :param tau_v: (default: 1)
        The optical depth at 5500\AA, used to normalize the
        attenuation curve.

    :returns tau:
        The optical depth at each wavelength.
    """
    if len(breaks) == len(alpha)+1:
        print("make sure of your power law breaks")
    tau = np.array(len(wave))
    for i in range(alpha):
        inds = (wave > breaks[i]) & (wave <= breaks[i+1])
        tau[inds] = tau_v * (wave / 5500)**alpha[i]
    return tau


def wg00(wave, tau_v=1, geometry='SHELL', composition='MW',
         local='homogenous', **kwargs):
    """ Witt+Gordon 2000 DIRTY radiative transfer results, for
    idealized geometries.
    """
    raise NotImplementedError

# ------------------
# EXTINCTION CURVES
# ------------------


def cardelli(wave, tau_v=1, R_v=3.1, **kwargs):
    r""" Cardelli, Clayton, and Mathis 1998 Milky Way extinction curve,
    with an update in the near-UV from O'Donnell 1994

    :param wave:
        The wavelengths at which optical depth estimates are desired.

    :param tau_v: (default: 1)
        The optical depth at 5500\AA, used to normalize the
        attenuation curve.

    :param R_v: (default: 3.1)
        The ratio of total selective extinction, parameterizing the
        slope of the attenuation curve.  A_v = R_v * E(B-V)

    :returns tau:
        The optical depth at each wavelength.
    """
    # if (wave < 1e3).any() :
    #     warnings.warn('Cardelli: extinction not defined (set to zero) below 1000AA')
    mic = wave*1e-4
    x_sup, x_inf = 10.0, 0.3
    x = 1 / mic
    a = np.zeros_like(x)
    b = np.zeros_like(x)

    w1 = (x >= 1.1) & (x <= 3.3)  # Optical 0.303 to 0.909 micron
    w2 = (x >= x_inf) & (x < 1.1)  # NIR 0.909 to 3.3 micron
    w3 = (x > 3.3) & (x <= 8)  # UV 0.125 - 0.303 micron
    w4 = (x > 8.0) & (x <= x_sup)  # XUV, 1000 -1250AA
    wsh = x > x_sup
    wlg = x < x_inf

    y = x[w1] - 1.82
    a[w1] = (1 + 0.17699 * y - 0.50447 * y**2. - 0.02427 * y**3. +
             0.72085 * y**4. + 0.01979 * y**5. - 0.77530 * y**6. +
             0.32999 * y**7.0)
    b[w1] = (1.41338 * y + 2.28305 * y**2. + 1.07233 * y**3. -
             5.38434 * y**4. - 0.62251 * y**5. + 5.30260 * y**6. -
             2.09002 * y**7.)

    y = x[w2]**1.61
    a[w2] = 0.574 * y
    b[w2] = -0.527 * y

    fa = x[w3] * 0.
    fb = x[w3] * 0.
    ou = (x[w3] > 5.9)
    # print(type(ou),ou[0], type(w3))

    if ou.any():
        y = x[w3][ou] - 5.9
        fa[ou] = -0.04473 * y**2. - 0.009779 * y**3.
        fb[ou] = 0.2130 * y**2. + 0.1207 * y**3.
    a[w3] = 1.752 - 0.316 * x[w3] - 0.104 / ((x[w3] - 4.67)**2. + 0.341) + fa
    b[w3] = -3.090 + 1.825 * x[w3] + 1.206 / ((x[w3] - 4.62)**2. + 0.263) + fb

    y = x[w4] - 8.
    a[w4] = -1.073 - 0.628 * y + 0.137 * y**2. - 0.070 * y**3.
    b[w4] = 13.670 + 4.257 * y - 0.420 * y**2. + 0.374 * y**3.

    tau = a + b / R_v
    return tau_v * tau


def smc(wave, tau_v=1, **kwargs):
    r"""Pei 1992 SMC extinction curve.

    :param wave:
        The wavelengths at which optical depth estimates are desired.

    :param tau_v: (default: 1)
        The optical depth at 5500\AA, used to normalize the
        attenuation curve.

    :returns tau:
        The optical depth at each wavelength.
    """
    if (wave < 1e3).any():
        warnings.warn('SMC: extinction extrapolation below 1000AA is poor')
    mic = wave * 1e-4
    aa = [185.,  27.,  0.005, 0.010, 0.012, 0.030]
    ll = [0.042, 0.08, 0.22,  9.7,   18.,   25.]
    bb = [90.,   5.50, -1.95, -1.95, -1.80, 0.00]
    nn = [2.0,   4.0,  2.0,   2.0,   2.0,   2.0]

    abs_ab = np.zeros_like(mic)
    norm_v = 0  # hack to go from tau_b to tau_v
    mic_5500 = 5500 * 1e-4

    for i, a in enumerate(aa):
        norm_v += aa[i] / ((mic_5500 / ll[i])**nn[i] +
                           (ll[i] / mic_5500)**nn[i] + bb[i])
        abs_ab += aa[i] / ((mic / ll[i])**nn[i] + (ll[i] / mic)**nn[i] + bb[i])

    return tau_v * (abs_ab / norm_v)


def lmc(wave, tau_v=1, **kwargs):
    r""" Pei 1992 LMC extinction curve.

    :param wave:
        The wavelengths at which optical depth estimates are desired.

    :param tau_v: (default: 1)
        The optical depth at 5500\AA, used to normalize the
        attenuation curve.

    :returns tau:
        The optical depth at each wavelength.
    """
    if (wave < 1e3).any():
        warnings.warn('LMC: extinction extrapolation below 1000AA is poor')
    mic = wave * 1e-4
    aa = [175.,  19.,  0.023, 0.005, 0.006, 0.020]
    ll = [0.046, 0.08, 0.22,  9.7,   18.,   25.]
    bb = [90.,   5.50, -1.95, -1.95, -1.80, 0.00]
    nn = [2.0,   4.5,  2.0,   2.0,   2.0,   2.0]

    abs_ab = mic * 0.
    norm_v = 0  # hack to go from tau_b to tau_v
    mic_5500 = 5500 * 1e-4

    for i, a in enumerate(aa):
        norm_v += aa[i] / ((mic_5500 / ll[i])**nn[i] +
                           (ll[i] / mic_5500)**nn[i] + bb[i])
        abs_ab += aa[i] / ((mic / ll[i])**nn[i] + (ll[i] / mic)**nn[i] + bb[i])

    return tau_v * (abs_ab / norm_v)
