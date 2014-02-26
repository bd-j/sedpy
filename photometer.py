# Simple Aperture photometry

import numpy as np
from numpy import hypot, sqrt
from astropy import wcs
import sys
import sklearn
from astroML.density_estimation import bayesian_blocks
from scipy.optimize import curve_fit

import matplotlib.pyplot as pl

thismod = sys.modules[__name__]

class Aperture(object):

    def world_to_pixels(shape, wcs):
        
        pass

    def object_flux(self, shape, image, ivar = None):
        """Measure total flux within an aperture (source + background)"""
        inds, fracs = self.pixnum(**shape)
        unc = 0
        if ivar is not None:
            unc = sqrt((fracs/ivar[inds[0], inds[1]]).sum())
        return (image[inds[0], inds[1]]*fracs).sum(), fracs.sum(), unc

    def measure_flux(self, shape, image, wcs = None, skypars = None, ivar = None):
        """Measure background subtracted flux."""
        o, a, e = self.object_flux(shape, image, ivar = ivar)
        b, ba, be = self.background.evaluate(image, skypars)
        flux = o - a*b
        flux_var = e*e + a*a*be*be/ba
        return flux, sqrt(flux_var)

    #def get_flux(self, image, ivar = None):
    #    return self.measure_flux(self.shape, image, ivar = ivar, skypars = self.skypars, wcs = self.wcs)
        
    
class Circular(Aperture):

    def __init__(self,  exact = False):
        if exact is True:
            self.pixnum = circle_frac_exact
        else:
            self.pixnum = circle_frac_quick
        self.background = ZeroSky()

class Elliptical(Aperture):

    def __init__(self):
        self.pixnum = ellipse_frac_quick()
        self.background = ZeroSky()

class Box(Aperture):
    
    def __init__(self):
        self.pixnum = box_frac_quick()
        self.background = ZeroSky()
#
######## Classes for sky measurement ######

class Background(object):

    def evaluate(self, image, skypars):
        inds, fracs = self.pixnum(**skypars)
        value, sdpp = self.skystats(image[inds[0], inds[1]], **skypars)
        return value, len(inds), sdpp

class Annulus(Background):

    def __init__(self, bgtype = 'quartile_sky'):
        self.pixnum = circle_frac_quick
        self.skystats = getattr(thismod, bgtype)

class EllipticalAnnulus(Background):

    def __init__(self, bgtype = 'quartile_sky'):
        self.pixnum = ellipse_frac_quick
        self.skystats = getattr(thismod,bgtype)

        
class ZeroSky(Background):
    """A class for sky values of zero, or for user defined sky statistics.
    The return_value is a tuple giving (sky, sky_area, sigma_sky_per_pixel)"""
    
    def __init__(self, bgtype = 'quartile_sky', return_value = (0,1,0)):
        self.pixnum = None
        self.skystats = None
        self.return_value = return_value

    def evaluate(self,image, skypars):
        return self.return_value


    
### Pixnum methods ####
def circle_frac_quick(xcen = 0, ycen = 0, radius = 1, inner_radius = None, subpixels = 1, **extras):
    """obtain fractional pixel coverage.  optionally use subpixels to
    increase precison (though this doesn't seem to help).  Assumes pixel
    centers have coordinates X.5, Y.5 """
    #setup
    center = np.array([xcen,ycen])
    sz = np.ceil((radius+1)*2)
    start = np.floor(center +0.5-radius)
    center = center*subpixels
    radius  = radius*subpixels
    sz = sz*subpixels
    start = (start-1)*subpixels
    if (start < 0).any():
        raise ValueError('Aperture extends off image edge')
    off = center - start - 0.5
    yy, xx = np.ogrid[ 0:sz, 0:sz ]
    rr = hypot(xx - off[0], yy-off[1])

    #find pixels within the radius
    within = (radius+0.5) - rr
    within[within > 1.0] = 1.0
    within[within < 0] = 0.

    #if it's an annulus
    if inner_radius is not None:
        within_inner = inner_radius*subpixels + 0.5 - rr
        within_inner[within_inner < 0.0] = 0.0
        within_inner[within_inner > 1.0] = 1.0
        within = within - within_inner
    an = within

    #rebin if you used subpixels
    if subpixels != 1:
        an = an.reshape((an.shape[0]/subpixels, subpixels, an.shape[1]/subpixels, subpixels)).mean(1).mean(2)
    #pick the pixels to rturn, and get their fractional coverage
    pix1 = np.where(an > 0.0)
    fracs = an[pix1[0],pix1[1]]
    x = (pix1[0] + start[0]/subpixels).astype('i8')
    y = (pix1[1] + start[1]/subpixels).astype('i8')

    return (x, y), fracs


def circle_frac_exact(xcen, ycen, radius):
    pass


def ellipse_frac_quick(xcen = 0, ycen = 0, a = 1, b = 1, pa = 90, precision = None):
    pass


#####SKY statistics determination methods #####


def quartile_sky(values, percentiles = [0.16, 0.5, 0.84], **extras):
    """Use the median and 16th percentile to estimate the standard
    deviation per pixel."""
    
    percentiles = np.asarray(percentiles)
    npix = len(values)
    #oo = np.argsort(values)
    qval = np.sort(values)[np.round(npix*percentiles).astype('i8')]
    #qval = values[oo[np.round(npix*percentiles)]] 
    return qval[1], qval[1]-qval[0]

def gaussfit_sky(values, p_thresh = 0.65, **extras):
    """Fit a gaussian to the lower part of a histogram of the sky values.
    The histogram bins are estimated using Bayesian blocks.  p_thresh gives
    the percentile below which the gaussian is fitted to the data. Return
    central value and estimate of standard deviation per pixel """
    
    bins = bayesian_blocks(values)
    print(len(bins),bins)
    #dbin = bins[1:]-bins[:-1]
    cbin = (bins[1:]+bins[:-1])/2
    hist = np.histogram(values, bins = bins, range = (bins.min(), bins.max()), density = True)
    
    #pdf = hist/dbin
    val_thresh = np.percentile(values, p_thresh)
    lower = cbin < p_thresh

    def gauss(x, *p):
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))

    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [np.max(hist[0]), values.mean(), values.std()]
    coeff, var_matrix = curve_fit(gauss, cbin[lower], hist[0][lower], p0=p0)
    print(len(hist[1]), len(hist[0]),type(coeff))
    pl.figure()
    pl.plot(cbin,hist[0], color = 'b')
    pl.plot(cbin, gauss(cbin, [coeff[0], coeff[1], coeff[2]]), color = 'r')
    pl.axvline(val_thresh)
    return coeff[1], coeff[2]

def gmm_sky(values, **extras):
    """Use a gaussian mixture model, via expectation maximization.
    of course, there's only one gaussian.  could add another for
    faint sources, bad pixels, but..."""
    
    gmm = sklearn.mixture.GMM()
    r = gmm.fit(values)
    return r.means_[0, 0], np.sqrt(r.covars_[0, 0])

def sigclip_sky(values, sigma = [3, 2.25], minlength = 5, **extras):
    """Use iterative sigma clipping"""
    
    def between(vals, sigs):
        m, s = vals.mean(), vals.std()
        return (vals < m+sig[1]*s) & (vals > m-sig[0]*s)
    
    while ( (False in between(values, sigma)) & (len(values) > minlength) ):
        values = values[between(values,sigma)]

    return values.mean(), values.std()


##### Centroiding #######
def centroid(images):
    """Dumb dumb centroiding.  assumes x and y axes are the
    last two dimensions of images.  Something is wrong with the
    broadcasting.  absolutely *have* to include weighting"""
    sz = images.shape[-2:]
    xg = np.arange(sz[0])
    yg = np.arange(sz[1])
    denom = images.sum(axis = (-1, -2))
    y = (yg[None,None,:]*images).sum(axis = (-2, -1)) / denom
    x = (xg[None,:,None]*images).sum(axis = (-2, -1)) / denom
    return x, y
