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

class FMcurve(object):
    """Class for handling extinction curves with the parameterization
    of Fitzpatrick & Massa 1990.  coeffs() method should be defined
    separately for each subclass."""
 
    def acurve(self, wave, A_v =1, R_v = None, f_bump = 1, **extras):
        x = 1e4/wave #inverse microns
        if R_v is not None:
            self.R_v = R_v
        if type(self.R_v) is not np.ndarray:
            self.R_v = np.atleast_1d(np.array(self.R_v))
        self.coeffs()
        return self.fm_curve(x, A_v, self.R_v, f_bump)

    def fm_curve(self, x, A_v, R_v, f_bump):
        return (A_v/R_v) * (self.c[0] + self.c[1]*x +
                       f_bump*self.c[2]*self._D(x, self.gamma, self.x0) +
                       self.c[3]*self._F(x) ) + A_v

    
    def _D(self, x, gamma, x0):
        """Drude profile for the 2175AA bump"""
        return x**2 / ( (x**2-x0**2)**2 +(x*gamma)**2 )

    def _F(self, x):
        """UV rise in the extinction curve"""
        return (x > 5.9)*(0.5392*(x-5.9)**2+0.05644*(x-5.9)**3)

    def coeffs(self):
        #c: linear intercept, linear slope, bump strength, UV rise 
        self.c = np.array([-0.0677, 0.698, 3.23, 0.400])
        self.c_unc = np.array([0.0, 0.0, 0.0, 0.0])
        self.x0 = 4.598  #bump center (inverse microns)
        self.x0_unc = 0.000
        self.gamma = 1.00 #bump width (inverse microns)
        self.gamma_unc = 0.000


class MilkyWay(FMcurve):
    """R dependent MW curves from Fitzpatrick 1999, including a variable bump (c_3)"""
    def __init__(self, R_v = 3.1):
        self.R_v = R_v
        self.coeffs()

    def coeffs(self):
        """FM99 parameterization coefficients."""
        if type(self.R_v) is not np.ndarray:
            self.R_v = np.atleast_1d(np.array(self.R_v))
        c2 = -0.824 + 4.717/self.R_v
        c1 = 2.030 - 3.007 * c2
        c3 = np.zeros(c1.shape)+3.23
        c4 = np.zeros(c1.shape)+0.41
        self.c = np.array([c1, c2, c3, c4])
        self.c_unc = np.array([0.,0.,0., 0.])
        self.x0 = np.zeros(c1.shape)+4.596
        self.x0_unc = 0. #not given
        self.gamma = np.zeros(c1.shape)+0.99
        self.gamma_unc = 0. #not given

    def acurve(self, wave, A_v =1, R_v = None, f_bump = 1, **extras):
        x = 1e4/wave #inverse microns
        if R_v is None:
            R_v = self.R_v
        else:
            self.R_v = R_v
            #print('Refreshing coeffs for R_v={0}'.format(R_v))
            self.coeffs()

        return np.where(x < 1e4/2700.,
                          self.optical(x,A_v,R_v, f_bump),
                          self.fm_curve(x, A_v, R_v, f_bump))
            

    def optical(self, x, A_v, R_v, f_bump):
        ebv = A_v/R_v
        
        # From Table 3 of Fitzpatrick 1999, for R_v = 3.1
        spline_x = 1e4/np.array([1,26500, 12200, 6000, 5470, 4670., 4100., 2700., 2600.])
        spline_x[0] = 0
        spline_a = np.array([0.,0.265,0.829, 2.688, 3.055, 3.806, 4.315, 6.265, 6.591])
        #Fix for differnt R_v values
        spline_a[0:3] *= R_v/3.1 #ir region
        # From Table 4
        spline_a[3:7] = np.array([-0.426+1.0044*R_v,
                                  -0.050+1.0016*R_v,
                                  0.701 +1.0016*R_v,
                                  1.208+1.0032*R_v-0.00033*R_v**2]) #optical
        spline_a[7:] = self.fm_curve(spline_x[7:], A_v, R_v, f_bump)/ebv #uv
        sp = interp.InterpolatedUnivariateSpline(spline_x, spline_a)
        return sp(x)*ebv
        
        


class LMC(FMcurve):
    """LMC average from Gordon et al. 2003. """
    def __init__(self, R_v = np.atleast_1d(np.array([3.41])) ):
        self.R_v = R_v
        self.coeffs()
    

    def coeffs(self):
        #c2 = -0.824 + 4.717/self.R_v
        #c1 = 2.030 - 3.007 * c2

        self.c = np.array([-0.890, 0.998, 2.719, 0.400])
        #self.c = np.array([c1, c2, 2.719, 0.400])
        self.c_unc = np.array([0.142, 0.027, 0.137, 0.036])
        self.x0 = 4.579
        self.x0_unc = 0.007
        self.gamma = 0.934
        self.gamma_unc = 0.016

class LMCshell(FMcurve):
    """LMC supershell average from Gordon et al. 2003"""
    def __init__(self, R_v = np.atleast_1d(np.ndarray([2.76])) ):
        self.R_v = R_v
        self.coeffs()

    def coeffs(self):
        self.c = np.array([-1.475, 1.132, 1.463, 0.294])
        self.c_unc = np.array([0.152, 0.029, 0.121, 0.057])
        self.x0 = 4.558
        self.x0_unc = 0.021
        self.gamma = 0.945
        self.gamma_unc = 0.026

class SMC(FMcurve):
    """SMC bar average from Gordon et al. 2003"""
    def __init__(self):
        self.coeffs()

    def coeffs(self):
        self.c = np.array([-4.959, 2.264, 0.389, 0.461])
        self.c_unc = np.array([0.197, 0.040, 0.110, 0.079])
        self.x0 = 4.60
        self.x0_unc = 0.0 #was held fixed
        self.gamma = 1.00
        self.gamma_unc = 0.0 #was held fixed


###### ATTENUATION CURVES ###########
######

def chevallard(wave, tau_v = 1, **kwargs):
    """\tau_v dependent attenuation curves matched to disk RT models, as in
    Chevallard et al. 2013.  No UV bump (or indeed tests in the UV at all)
    """
    #missing a UV bump
    alpha_v = 2.8/(1+np.sqrt(tau_v)) #+/- 25%
    bb = 0.3 - 0.05*tau_v  #+/- 10%
    alpha = alpha_v + bb*(wave*1e-4 - 0.55)
    tau_lambda = tau_v*(wave/5500.0)**(0-alpha)
    return tau_lambda

def powerlaw(wave, tau_v = 1, alpha = 1.0, **kwargs):
    return tau_v * (wave/5500)**(0-alpha)

def calzetti(wave, tau_v = 1, R_v=4.05, **kwargs):
    """Calzetti et al. 2000 starburst attenuation curve, with extrapolations
    to the FUV and NIR"""
    p11=1/0.11
    ff11=2.659*(-2.156+1.509*p11-0.198*p11**2.+0.011*p11**3.0) + R_v
    p12=1/0.12
    ff12=2.659*(-2.156 + 1.509*p12 - 0.198*p12**2. + 0.011*p12**3) + R_v
    slope1=(ff12 - ff11)/100.
    ff99=2.659*(-1.857 + 1.040/2.19) + R_v
    ff100=2.659*(-1.857 + 1.040/2.2) + R_v
    slope2=(ff100 - ff99)/100.
    #do it
    x=1e4/wave
    ff = ( (wave >= 6300.) & (wave <= 22000) ) * (2.659*(-1.857 + 1.040*x) + R_v)
    ff +=( ( (wave >= 1200.) & (wave < 6300) ) *
           (2.659*(-2.156 + 1.509*x - 0.198*x**2. + 0.011*x**3.) + R_v) )
    ff += (wave < 1200.) * (ff11 + (wave-1100.)*slope1)
    ff += (wave > 22000.) * (ff99 + (wave-21900.)*slope2)

    ff[ff < 0] = 0
    tau_lambda=tau_v*ff/R_v/0.999479
    return tau_lambda


def wg00(wave,tau_v =1, geometry = 'SHELL',composition='MW', local = 'homogenous', **kwargs):
    """Witt+Gordon 2000 DIRTY radiative transfer results, for idealized geometries"""
    pass

def conroy(wave, tau_v =1, R_v=3.1, f_bump=0.6, **kwargs):
    """Conroy et al 2010 dust attenuation curves including a decreased UV bump."""
    x = 1e4/wave
    nx = x.shape[0]
    a = x * 0
    b = x * 0 

    #IR 0.909 - 3.3 micron
    ir = np.where( (x >= 0.3) & (x < 1.1) )
    a[ir] = 0.574*(x[ir]**1.61)
    b[ir] = (-0.527)*(x[ir]**1.61)

    #optical 0.303 - 0.909 micron
    opt = np.where( (x >= 1.1) & (x < 3.3) )
    y = x[opt]-1.82
    a[opt] = ( 1 + 0.177*y - 0.504*y**2-0.0243*y**3 + 0.721*y**4 +
	      0.0198*y**5 - 0.775*y**6 + 0.330*y**7 )
    b[opt] = ( 1.413*y + 2.283*y**2 + 1.072*y**3 - 5.384*y**4 -
	       0.622*y**5 + 5.303*y**6 - 2.090*y**7 )

    #NUV 0.17 to 0.303 micron
    nuv = np.where( (x >= 3.3) & (x < 5.9) )
    fa = (3.3/x[nuv])**6. * (-0.0370+0.0469*f_bump-0.601*f_bump/R_v+0.542/R_v)
    a[nuv] = 1.752 - 0.316*x[nuv] - (0.104*f_bump/((x[nuv]-4.67)**2+0.341)) + fa
    b[nuv] = (-3.09) + 1.825*x[nuv] + (1.206*f_bump/((x[nuv]-4.62)**2+0.263))

    #FUV 0.125 - 0.17 micron
    fuv = np.where( (x >= 5.9) & (x < 8.0) )
    fa = (-0.0447)*(x[fuv]- 5.9)**2.0 - 0.00978*(x[fuv]- 5.9)**3.
    fb = 0.213*(x[fuv]- 5.9)**2. + 0.121*(x[fuv]- 5.9)**3.
    a[fuv] = 1.752 - 0.316*x[fuv] - (0.104*f_bump/((x[fuv]-4.67)**2+0.341)) + fa
    b[fuv] = (-3.09) + 1.825*x[fuv] + (1.206*f_bump/((x[fuv]-4.62)**2+0.263)) + fb

    alam = (a+b/R_v)

    #XUV below 1250AA
    xuv = np.where(x >= 8.0)
    x8=8.0
    fa = (-0.0447)*(x8- 5.9)**2.-0.00978*(x8- 5.9)**3.
    fb = 0.213*(x8- 5.9)**2. + 0.121*(x8- 5.9)**3.
    af = 1.752 - 0.316*x8 - (0.104*f_bump/((x8-4.67)**2+0.341)) + fa
    bf = (-3.09) + 1.825*x8 + (1.206*f_bump/((x8-4.62)**2+0.263)) + fb
    a8 = (af+bf/R_v)

    alam[xuv] = (x8/x[xuv])**(-1.3)*a8

    return tau_v*alam

def wild_powerlaw(wave, tau_v =1 ,alpha=[0.7,0.7,0.7], breaks=[0,3000,10000,4e4], **kwargs):
    """As in V. Wild 2011, i.e. power-law slope can change between regions.
    Superceded by Chevallard 2013 for optical/NIR"""
    if len(breaks) == len(alpha)+1 : print "make sure of your power law breaks"
    tau=np.array(len(wave))
    for i in range(alpha):
	inds=np.where((wave > breaks[i]) & (wave <=breaks[i+1]))
	tau[inds]=tau_v*(wave/5500)**alpha[i]
    return tau


##### EXTINCTION CURVES ########

def cardelli(wave, tau_v = 1, R_v=3.1, **kwargs):
    """Cardelli, Clayton, and Mathis 1998 Milky Way extinction curve,
    with an update in the near-UV from O'Donnell 1994"""

    if (wave < 1e3).any() :
        warnings.warn('Cardelli: extinction not defined (set to zero) below 1000AA')
    mic=wave*1e-4
    x_sup, x_inf = 10.0, 0.3
    x = 1/mic
    a = x * 0.
    b = x * 0.

    w1 = np.where( (x >= 1.1) & (x <= 3.3) ) #Optical 0.303 to 0.909 micron
    w2 = np.where( (x >= x_inf) & (x < 1.1) ) # NIR 0.909 to 3.3 micron
    w3 = np.where( (x > 3.3) & (x <= 8) ) #UV 0.125 - 0.303 micron
    w4 = np.where( (x > 8.0) & (x <= x_sup) ) #XUV, 1000 -1250AA
    wsh = np.where( x > x_sup)
    wlg = np.where( x < x_inf)

    print(w2[0])

    y = x[w1] - 1.82
    a[w1] = 1. + 0.17699 * y - 0.50447 * y**2. - 0.02427 * y**3. + 0.72085 * y**4. + \
	0.01979 * y**5. - 0.77530 * y**6. + 0.32999 * y**7.0
    b[w1] = 1.41338 * y + 2.28305 * y**2. + 1.07233 * y**3. - 5.38434 * y**4. - \
	0.62251 * y**5. + 5.30260 * y**6. - 2.09002 * y**7.

    y = (x[w2])**1.61
    a[w2] = 0.574 * y
    b[w2] = -0.527 * y

    fa = x[w3] * 0.
    fb = x[w3] * 0.
    ou = (x[w3] > 5.9)
    print(type(ou),ou[0], type(w3))
    
    if ou.any():
	y = x[w3[0][ou]] - 5.9
	fa[ou] = -0.04473 * y**2. - 0.009779 * y**3.
	fb[ou] = 0.2130 * y**2. + 0.1207 * y**3.
    a[w3] = 1.752 - 0.316 * x[w3] - 0.104 / ((x[w3] - 4.67)**2. + 0.341) + fa
    b[w3] = -3.090 + 1.825 * x[w3] + 1.206 / ((x[w3] - 4.62)**2. + 0.263) + fb

    y = x[w4] - 8.
    a[w4] = -1.073 - 0.628 * y + 0.137 * y**2. - 0.070 * y**3.
    b[w4] = 13.670 + 4.257 * y - 0.420 * y**2. + 0.374 * y**3.

    tau = a + b / R_v
    return tau_v*tau

def smc(wave, tau_v =1, **kwargs):
    """Pei 1992 SMC extinction curve"""
    if (wave < 1e3).any() :
        warnings.warn('SMC: extinction extrapolation below 1000AA is poor')
    mic = wave * 1e-4
    aa = [185.,  27.,  0.005, 0.010, 0.012, 0.030]
    ll = [0.042, 0.08, 0.22,  9.7,   18.,   25.]
    bb = [90.,   5.50, -1.95, -1.95, -1.80, 0.00]
    nn = [2.0,   4.0,  2.0,   2.0,   2.0,   2.0]

    abs_ab = np.zeros_like(mic)
    norm_v = 0 #hack to go from tau_b to tau_v
    mic_5500=5500*1e-4

    for i, a in enumerate(aa):
	norm_v +=  aa[i] / ( (mic_5500/ll[i])**nn[i] + (ll[i]/mic_5500)**nn[i] + bb[i] ) 
	abs_ab +=  aa[i] / ((mic/ll[i])**nn[i] + (ll[i]/mic)**nn[i] + bb[i])

    return tau_v * (abs_ab/norm_v)

def lmc(wave, tau_v =1, **kwargs):
    """Pei 1992 LMC extinction curve"""
    if (wave < 1e3).any() :
        warnings.warn('LMC: extinction extrapolation below 1000AA is poor')
    mic = wave * 1e-4
    aa = [175.,  19.,  0.023, 0.005, 0.006, 0.020]
    ll = [0.046, 0.08, 0.22,  9.7,   18.,   25.]
    bb = [90.,   5.50, -1.95, -1.95, -1.80, 0.00]
    nn = [2.0,   4.5,  2.0,   2.0,   2.0,   2.0]

    abs_ab = mic * 0.
    norm_v = 0 #hack to go from tau_b to tau_v
    mic_5500=5500*1e-4

    for i, a in enumerate(aa):
	norm_v +=  aa[i] / ( (mic_5500/ll[i])**nn[i] + (ll[i]/mic_5500)**nn[i] + bb[i] ) 
	abs_ab +=  aa[i] / ((mic/ll[i])**nn[i] + (ll[i]/mic)**nn[i] + bb[i])

    return tau_v * (abs_ab/norm_v)

    

