# Python module for storing filter information and tools for
# projecting spectra onto filters.  Also includes tools for convolving
# spectra.
#
# Assumed input units are erg/s/cm^2/AA and AA

#To Do: -- Add methods for obtaining magnitude uncertainties
#          given an error spectrum along with the source spectrum
#       -- Add index measurements
#       -- Test the broadening functions for accuracy and speed.
#       -- Add some redshifting methods/classes?


import numpy as np
import os
try:
    from pkg_resources import resource_filename
except ImportError:
    pass
from yanny import read as yanny_read
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

lightspeed = 2.998e18 #AA/s
## Load useful reference spectra ######
sedpydir, f = os.path.split(__file__)
sedpydir = sedpydir
try:
    vega_file = resource_filename('data', 'alpha_lyr_stis_005.fits')
except:
    vega_file = os.path.join(sedpydir, 'data','alpha_lyr_stis_005.fits')
    
#this file should be in AA and erg/s/cm^2/AA
if os.path.isfile( vega_file ):
    fits = pyfits.open( vega_file )
    vega = np.column_stack( (fits[1].data.field('WAVELENGTH'), fits[1].data.field('FLUX')) )
    fits.close()
else:
    raise ValueError('Could not find Vega spectrum at {0}'.format(vega_file))

try:
    solar_file = resource_filename('data','sun_kurucz93.fits')
except:
    solar_file = os.path.join(sedpydir,'data','sun_kurucz93.fits')

rat = (1.0/(3600*180/np.pi*10))**2.0 # conversion to d=10 pc from 1 AU
#this file should be in AA and erg/s/cm^2/AA at 1AU
if os.path.isfile( solar_file ):
    fits = pyfits.open( solar_file )
    solar = np.column_stack( (fits[1].data.field('WAVELENGTH'), fits[1].data.field('FLUX')*rat) )
    fits.close()
else:
    raise ValueError('Could not find Solar spectrum at {0}'.format(solar_file))


class Filter(object):
    """
    This class operates on filter transmission files.  It reads
    SDSS-style yanny files containing filter transmissions (these are
    easy to create) and determines a number of useful filter
    quantities.  Methods are provided to convolve a source spectrum
    with the filter and return the magnitude.

    :param kname: (default: 'sdss_r0')
        The kcorrect style name of the filter, excluing '.par',
        e.g. sdss_r0.

    :param nick: (optional)
        A nickname to associate with the filter.
        
    """

    ab_gnu=3.631e-20   #AB reference spctrum in erg/s/cm^2/Hz
    npts=0
    
    def __init__(self, kname = 'sdss_r0', nick = None):
        """
        Constructor.
        """
        self.name = kname
        if nick is None :
            self.nick = kname
        else:
            self.nick = nick

        try:
            self.filename = resource_filename('/data/filters/',kname + '.par')
        except:
            self.filename = sedpydir + '/data/filters/' + kname + '.par'
        if type( self.filename ) == type( '' ):
            if not os.path.isfile( self.filename ):
                raise ValueError( 'Filter transmission file {0} does not exist!'.format(self.filename) )
            self.loadKFilter(self.filename)

    def loadKFilter(self, filename):
        """
        Read a filter in kcorrect (yanny) format and populate the
        wavelength and transmission arrays.  Then determine a number
        of filter properties and store in the object.

        :param filename:
            The fully qualified path and filename of the yanny file
            that contains the filter transmission.
        """

        ##This should be replaced with the sdsspy yanny file readers
        ##Done, code kept here in case reversion required
        #f=open(filename,'rU')
        #wave=[]
        #trans=[]
        #for line in f:
        #    cols=line.split()
        #    if len(cols) > 2 and cols[0].find('KFILTER') > -1:
        #        wave.append(float(cols[1]))
        #        if cols[0].find('SDSS') > -1:
        #            trans.append(float(cols[4])) #use airmass=1.3 passband.  HACKY
        #        else:
        #            trans.append(float(cols[2]))
        #f.close()

        ff = yanny_read(filename, one = True)
        wave = ff['lambda']
        trans = ff['pass']
        #clean negatives, NaNs, and Infs, then sort, then store
        ind = np.where(np.logical_and( np.isfinite(trans), (trans >= 0.0) ))[0]
        order = wave[ind].argsort()
        self.npts = ind.shape[0]
        self.wavelength = wave[ind[order]]
        self.transmission = trans[ind[order]]

        self.getProperties()

        
    def getProperties(self):
        """
        Determine and store a number of properties of the filter and
        store them in the object.  These properties include several
        'effective' wavelength definitions and several width
        definitions, as well as the in-band absolute AB solar
        magnitude, the Vega and AB reference detector signal, and the
        conversion between AB and Vega magnitudes.
        """

        i0 = np.trapz(self.transmission*np.log(self.wavelength), np.log(self.wavelength))
        i1 = np.trapz(self.transmission, np.log(self.wavelength))
        i2 = np.trapz(self.transmission*self.wavelength, self.wavelength)
        i3 = np.trapz(self.transmission, self.wavelength)
        i4 = np.trapz(self.transmission*( np.log(self.wavelength) )**2.0, np.log(self.wavelength))
        
        self.wave_effective = np.exp(i0/i1)
        self.wave_pivot = np.sqrt(i2/i1)
        self.wave_mean = self.wave_effective
        self.wave_average = i2/i3

        self.gauss_width = (i4/i1)**(0.5)
        self.effective_width = 2.0*np.sqrt( 2.0*np.log(2.0) )*self.gauss_width*self.wave_mean
        self.rectangular_width = i3/self.transmission.max()
        #self.norm  = np.trapz(transmission,wavelength)

        self.ab_counts = self.objCounts( self.wavelength,self.ab_gnu*lightspeed/(self.wavelength**2) )
        self.vega_counts = self.objCounts(vega[:,0],vega[:,1]) 
        self.ab_to_vega = -2.5*np.log10(self.ab_counts/self.vega_counts)
        if self. wave_mean < 1e5:
            self.solar_ab_mag = self.ABMag(solar[:,0],solar[:,1])
        else:
            self.solar_ab_mag = float('NaN')
            
    def display(self):
        """
        Plot the filter transmission curve.
        """
        if self.npts > 0:
            plt.plot(self.wavelength,self.transmission)
            plt.title(self.name)

    def objCounts(self, sourcewave, sourceflux, sourceflux_unc=0):
        """
        Project source spectrum onto filter and return the detector
        signal.

        :param sourcewave:
            Apectrum wavelength (in AA), ndarray of shape (nwave).
            
        :param sourceflux:
            Associated flux (assumed to be in erg/s/cm^2/AA), ndarray
            of shape (nspec,nwave).

        :returns counts:
            Detector signal(s) (nspec).
        """
        
        #interpolate filter transmission to source spectrum
        newtrans = np.interp(sourcewave, self.wavelength,self.transmission, 
                                left=0.,right=0.)
            #print(self.name,sourcewave.shape,self.wavelength.shape,newtrans.shape)
        
        #integrate lambda*f_lambda*R
        if True in (newtrans > 0.):
            ind = np.where(newtrans > 0.)
            ind=ind[0]
            counts = np.trapz(sourcewave[ind]*newtrans[ind]*sourceflux[...,ind], sourcewave[ind],axis=-1)
            #if  np.isinf(counts).any() : print(self.name, "Warn for inf value")
            return np.squeeze(counts)
        else:
            return float('NaN')

    def ABMag(self, sourcewave, sourceflux, sourceflux_unc=0):
        """
        Project source spectrum onto filter and return the AB
        magnitude,
        
        :param sourcewave:
            Spectrum wavelength (in AA), ndarray of shape (nwave).
            
        :param sourceflux:
            Associated flux (assumed to be in erg/s/cm^2/AA), ndarray
            of shape (nobj,nwave).

        :returns mag:
            AB magnitude of the source.
        """
        
        return 0-2.5*np.log10(self.objCounts(sourcewave, sourceflux)/self.ab_counts)

    def vegaMag(self, sourcewave, sourceflux, sourceflux_unc=0):
        """
        Project source spectrum onto filter and return the Vega
        magnitude.
               
        :param sourcewave:
            Spectrum wavelength (in AA), ndarray of shape (nwave).
            
        :param sourceflux:
            Associated flux (assumed to be in erg/s/cm^2/AA), ndarray
            of shape (nobj,nwave).

        :returns mag:
            Vega magnitude of the source.
        """

        return 0-2.5*np.log10(self.objCounts(sourcewave, sourceflux)/self.vega_counts)        


###Useful utilities#####

def load_filters(filternamelist):
    """
    Given a list of filter names, this method returns a list of Filter
    objects.

    :param filternamelist:
        List of strings giving names of the filters.

    :returns filterlist:
        A list of filter objects.
    """
    
    filterlist = []
    for f in filternamelist:
        #print(f)
        filterlist.append(Filter(f))
    return filterlist

def getSED(sourcewave, sourceflux, filterlist):
    """
    Takes wavelength vector, a flux array and list of Filter objects
    and returns the SED in AB magnitudes.

    :param sourcewave:
        Spectrum wavelength (in AA), ndarray of shape (nwave).
        
    :param sourceflux:
        Associated flux (assumed to be in erg/s/cm^2/AA), ndarray of
        shape (nsource,nwave).
        
    :param filterlist:
        List of filter objects, of length nfilt.

    :returns sed:
        array of broadband magnitudes, of shape (nsource, nfilter).
    """

    if filterlist is None:
        return None
    sourceflux = np.atleast_2d(sourceflux)
    sedshape = [sourceflux.shape[0], len(filterlist)]
    sed = np.zeros(sedshape)
    for i,f in enumerate(filterlist):
        sed[:,i] = f.ABMag(sourcewave,sourceflux)
    return np.squeeze(sed)

def filter_dict(filterlist):
    fdict = {}
    for i,f in enumerate(filterlist):
        fdict[f.nick] = i
    return fdict

###Routines for spectra######

def Lbol(wave,spec,wave_min=90,wave_max = 1e6):
    """
    Calculate the bolometric luminosity of a spectrum or spectra.

    :param wave:
       The wavelength vector of length nwave.
       
    :param spec:
       The spectra, of shape (...,nsource, nwave).
       
    :param wave_min:
       Minimum wavelength for the integral.
       
    :param max_wave:
       Maximum wavelength for the integral

    :returns lbol:
       The bolometric luminosity, integrated from wave_min to
       wave_max.  Array of length (...nsource)
    """

    inds = np.where(np.logical_and(wave < wave_max, wave >= wave_min))
    return np.trapz(spec[...,inds[0]],wave[inds])

def air2vac(air):
    """
    Convert from in-air wavelengths to vacuum wavelengths.
    Based on Allen's Astrophysical Quantities.

    :param air:
        The in-air wavelengths.
        
    :returns vac:
        The corresponding vacuum wavelengths.
    """
    ss = 1e4/air
    vac = air * (1 + 6.4328e-5 + 2.94981e-2 / (146 - ss**2) + 2.5540e-4 / (41 - ss**2))
    return vac

def vac2air(vac):
    """
    Convert from vacuum wavelengths to in-air wavelengths.  Follows
    the SDSS statement of the IAU standard from Morton 1991 ApJS.

    vac2air(air2vac(wave)) yields wave to within 1 part in a million
    over the optical range.

    :param vac:
        The vacuum wavelengths.
        
    :returns vac:
        The corresponding in-air wavelengths.

    """
    
    air = vac / (1.0 + 2.735182e-4 + 131.4182/vac**2 + 2.76249e8 / vac**4)
    return air



def vel_broaden(sourcewave, sourceflux, sigma_in, sigma0=0,
                outwave = None, nsig = 5.0, minusewave = 0, maxusewave = 1e8):
    """
    Vectorized version of velocity broadening.  This can become very
    slow when memory constraints are reached (i.e. when nw_in * nw_out
    * nsource is large).  This should be rewritten to work with (fft)
    convolutions for speed, though that requires regular (log)
    velocity scale.
    """
    sourceflux = np.atleast_2d(sourceflux)
    #sigma after accounting for intrinsic sigma of the library.
    sigma = np.sqrt(sigma_in**2-sigma0**2)
    
    if outwave is None:
        outwave = sourcewave
    # Set up arrays and limits
    minw, maxw = outwave.min(), outwave.max()
    maxw  *= (1 + nsig * sigma/lightspeed)
    minw  *= (1 - nsig * sigma/lightspeed)
    use = (sourcewave > np.max([minw, minusewave])) & (sourcewave < np.min([maxw, maxusewave]))
    
    K = 1 / (sigma * np.sqrt(2 * np.pi)) #don't need since renormalizing weights
    wr = outwave[:,None] / sourcewave[None, use ]
    v = lightspeed/1e13 * (1 - wr)
    ee = np.exp(-0.5 * v**2 / sigma**2)
    #renormalize - why?
    ee /= np.trapz(ee, x = v, axis=-1)[:,None]
    #now, integrate over the velocities in the sourcewave direction
    flux = np.trapz( sourceflux[:,None,use] * ee[None,:,:], x= v[None,:,:], axis = -1)

    #H = lightspeed/1e13/sigma
    #flux = np.trapz( sourceflux[:,None,use] * np.exp(-0.5 * (H * (1 - outwave[:,None] / sourcewave[None, use ]))**2),
    #                 x= lightspeed/1e13 * (1 - outwave[None,:,None]/sourcewave[None,None,:]), axis = -1)
    
    return  flux #* K


def vel_broaden_fast(sourcewave, sourceflux, sigma):
    pass

    
def wave_broaden(sourcewave, sourceflux, fwhm, fwhm0 = 0, outwave = None,
                 nsig = 5.0, minusewave = 0, maxusewave = 1e8):
    """
    Vectorized version of wavelength broadening.  This can become very
    slow when memory constraints are reached (i.e. when nw_in * nw_out
    * nsource is large).  This should be rewritten to work with (fft)
    convolutions for speed.
    """
    sourceflux = np.atleast_2d(sourceflux)
    
    sigma = np.sqrt(fwhm**2. - fwhm0**2.) / 2.3548
    if outwave is None:
        outwave = sourcewave
    # Set up arrays and limits
    minw, maxw = outwave.min(), outwave.max()
    maxw  *= (1 + nsig * sigma)
    minw  *= (1 - nsig * sigma)
    use = (sourcewave > np.max([minw, minusewave])) & (sourcewave < np.min([maxw, maxusewave]))

    dl = outwave[:,None] - sourcewave[None, use]
    ee = np.exp(-0.5 * dl**2 / sigma**2)
    ee /= np.trapz(ee, x = dl, axis=-1)[:,None]
    flux = np.trapz( sourceflux[:,None,use] * ee[None,:,:], x= dl[None,:,:], axis = -1)

    return flux


def broaden(sourcewave, sourceflux, width, width0 = 0, stype = 'vel', **kwargs):
    if stype is 'vel':
        return vel_broaden(sourcewave, sourceflux, width, sigma0 = width0, **kwargs)
    elif stype is 'wave':
        return wave_broaden(sourcewave, sourceflux, width, fwhm0 = width0, **kwargs)


#    K=(sigma*sqrt(2.*!PI))^(-1.)
#    for iw in range(len(outwave)):
#        dl=(outwave[iw]-wave)**2.
#        this=np.where(dl < length*sigma**2.,count)
#        if count > 0:
#            ee=exp(-0.5*dl[this]/sigma^2)
#            broadflux[iw]=int_tabulated(wave[this],flux[this]*ee)
#    return,broadflux

#def redshift(sourcewave,sourceflux, z):

def selftest():
    """Compare to the values obtained from the K-correct code
    (which uses a slightly different Vega spectrum)"""

    filternames=['galex_FUV','sdss_u0','sdss_g0','sdss_r0','sdss_i0','spitzer_irac_ch2']
    weff_kcorr=[1528.0,3546.0,4669.6,6156.25,7471.57,44826.]
    msun_kcorr=[18.8462,6.38989,5.12388,4.64505,4.53257,6.56205]
    ab2vega_kcorr=[2.3457,0.932765,-0.0857,0.155485,0.369598,3.2687]

    filterlist=loadFilters(filternames)
    for i in range(len(filterlist)):
        print(filterlist[i].wave_effective, filterlist[i].solar_ab_mag, filterlist[i].ab_to_vega)
        assert abs(filterlist[i].wave_effective-weff_kcorr[i]) < weff_kcorr[i]*0.01
        assert abs(filterlist[i].solar_ab_mag-msun_kcorr[i]) < 0.05
        #assert abs(filterlist[i].ab_to_vega+ab2vega_kcorr[i]) < 0.05 #this fails because of the vega spectrum used by k_correct

def selftest_broaden():
    import fsps
    sps = fsps.StellarPopulation()
    ages = [1., 2., 3., 4., 5., 6., 7.]
    allspec = np.zeros( [len(ages), len(sps.wavelengths)])
    for i,tage in enumerate(ages):
        wave, spec = sps.get_spectrum(peraa = True, tage =tage)
        allspec[i,:] = spec

    outwave = wave[(wave > 1e3) & (wave < 1e4)]
    vbflux = vel_broaden( wave, allspec, 150., outwave = outwave )
    wbflux = wave_broaden( wave, allspec, 5., outwave = outwave )

        
