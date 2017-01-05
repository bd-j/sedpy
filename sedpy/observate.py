# Python module for storing filter information and tools for projecting spectra
# onto filters.  Also includes tools for convolving spectra.
#
# Assumed input units are erg/s/cm^2/AA and AA

import numpy as np
import os
try:
    from pkg_resources import resource_filename, resource_listdir
except ImportError:
    pass
from .yanny import read as yanny_read
from .reference_spectra import vega, solar, sedpydir

__all__ = ["Filter", "load_filters", "getSED",
           "air2vac", "vac2air", "Lbol"]

lightspeed = 2.998e18  # AA/s


class Filter(object):
    """This class operates on filter transmission files.  It reads SDSS-style
    yanny files (these are easy to create) or simple 2-column ascii files
    containing filter transmission curves. It caches a number of useful filter
    quantities.  Methods are provided to project a source spectrum onto the
    filter and return the magnitude, AB or Vega.

    :param kname: (default: 'sdss_r0')
        The kcorrect style name of the filter, excluing '.par',
        e.g. sdss_r0.

    :param nick: (optional)
        A nickname to associate with the filter.

    :param directory: (optional)
        The path to the directory containing the filter file.  If not given
        then the sedpy/data/filters/ directory will be searched.
    """
    ab_gnu = 3.631e-20  # AB reference spctrum in erg/s/cm^2/Hz
    npts = 0

    def __init__(self, kname='sdss_r0', nick=None, directory=None, **extras):
        """Constructor.
        """
        self.name = kname
        if nick is None:
            self.nick = kname
        else:
            self.nick = nick

        if directory is None:
            try:
                self.filename = resource_filename('sedpy', '/data/filters/' +
                                              kname + '.par')
            except:
                self.filename = os.path.join(sedpydir, '/data/filters/', kname + '.par')
        else:
            self.filename = os.path.join(directory, kname+'.par')

        if isinstance(self.filename, str):
            if not os.path.isfile(self.filename):
                raise ValueError('Filter transmission file {0} does '
                                 'not exist!'.format(self.filename))
            try:
                self.load_kfilter(self.filename)
            except:
                self.load_filter(self.filename)

        self.get_properties()

    def __repr__(self):
        return '{}({})'.format(self.__class__, self.name)
                        
    def load_kfilter(self, filename):
        """Read a filter in kcorrect (yanny) format and populate the wavelength
        and transmission arrays.

        :param filename:
            The fully qualified path and filename of the yanny file that
            contains the filter transmission.
        """

        ff = yanny_read(filename, one=True)
        wave = ff['lambda']
        trans = ff['pass']
        # Clean negatives, NaNs, and Infs, then sort, then store
        ind = np.where(np.isfinite(trans) & (trans >= 0.0))[0]
        order = wave[ind].argsort()
        self.npts = ind.shape[0]
        self.wavelength = wave[ind[order]]
        self.transmission = trans[ind[order]]

    def load_filter(self, filename):
        """Read a filter in simple two column ascii format and populate the
        wavelength and transmission arrays.  The first column is wavelength in
        AA and the second column is transmission (detector signal per photon)
        
        :param filename:
            The fully qualified path and filename of the file that contains the
            filter transmission.
        """
        wave, trans = np.genfromtxt(filename, usecols=(0,1), unpack=True)
        # Clean negatives, NaNs, and Infs, then sort, then store
        ind = np.isfinite(trans) & (trans >= 0.0)
        order = wave[ind].argsort()
        self.npts = ind.sum()
        self.wavelength = wave[ind][order]
        self.transmission = trans[ind][order]

    def get_properties(self):
        """Determine and store a number of properties of the filter and store
        them in the object.  These properties include several 'effective'
        wavelength definitions and several width definitions, as well as the
        in-band absolute AB solar magnitude, the Vega and AB reference
        zero-point detector signal, and the conversion between AB and Vega
        magnitudes.

        See Fukugita et al. (1996) AJ 111, 1748 for discussion and definition
        of many of these quantities.
        """
        # Calculate some useful integrals
        i0 = np.trapz(self.transmission * np.log(self.wavelength),
                      np.log(self.wavelength))
        i1 = np.trapz(self.transmission,
                      np.log(self.wavelength))
        i2 = np.trapz(self.transmission * self.wavelength,
                      self.wavelength)
        i3 = np.trapz(self.transmission,
                      self.wavelength)

        self.wave_effective = np.exp(i0 / i1)
        self.wave_pivot = np.sqrt(i2 / i1)
        self.wave_mean = self.wave_effective
        self.wave_average = i2 / i3
        self.rectangular_width = i3 / self.transmission.max()

        i4 = np.trapz(self.transmission *
                      (np.log(self.wavelength / self.wave_effective))**2.0,
                      np.log(self.wavelength))
        self.gauss_width = (i4 / i1)**(0.5)
        self.effective_width = (2.0 * np.sqrt(2. * np.log(2.)) *
                                self.gauss_width *
                                self.wave_effective)
        # self.norm  = np.trapz(transmission,wavelength)

        # Get zero points and AB to Vega conversion
        self.ab_zero_counts = self.obj_counts(self.wavelength,
                                              self.ab_gnu * lightspeed /
                                              self.wavelength**2)
        # If blue enough get AB mag of vega
        if self.wave_mean < 1e6:
            self.vega_zero_counts = self.obj_counts(vega[:,0], vega[:,1])
            self._ab_to_vega = -2.5 * np.log10(self.ab_zero_counts /
                                               self.vega_zero_counts)
        else:
            self.vega_zero_counts = float('NaN')
            self._ab_to_vega = float('NaN')
        # If blue enough get absolute solar magnitude
        if self.wave_mean < 1e5:
            self.solar_ab_mag = self.ab_mag(solar[:,0], solar[:,1])
        else:
            self.solar_ab_mag = float('NaN')

    @property
    def ab_to_vega(self):
        """The conversion from AB to Vega systems for this filter.  It has the
        sense

        :math:`m_{Vega} = m_{AB} + Filter().ab_to_vega`
        """
        return self._ab_to_vega
            
    def display(self, normalize=False, ax=None):
        """Plot the filter transmission curve.
        """
        if self.npts > 0:
            if ax is None:
                import matplotlib.pyplot as pl
                fig, ax = pl.subplots()
                ax.title(self.nick)
                fig.show()
            if normalize:
                ax.plot(self.wavelength, self.transmission / self.transmission.max())
            else:
                ax.plot(self.wavelength, self.transmission)
            return ax

    def obj_counts(self, sourcewave, sourceflux, sourceflux_unc=0):
        """Project source spectrum onto filter and return the detector signal.

        :param sourcewave:
            Spectrum wavelength (in AA), ndarray of shape (nwave).  Must be
            monotonic increasing.

        :param sourceflux:
            Associated flux (assumed to be in erg/s/cm^2/AA), ndarray of shape
            (nspec,nwave).

        :returns counts:
            Detector signal(s) (nspec).
        """
        assert sourcewave[1] > sourcewave[0], "``sourcewave`` not in ascending order."
        # Interpolate filter transmission to source spectrum
        newtrans = np.interp(sourcewave, self.wavelength, self.transmission,
                             left=0., right=0.)

        # Integrate lambda*f_lambda*R
        if True in (newtrans > 0.):
            positive = np.where(newtrans > 0.)[0]
            edge = (positive.min() < 1) | (positive.max() >= (len(newtrans)-1))
            # assert ~edge, "Source spectrum does not span filter."
            ind = slice(max(positive.min() - 1, 0),
                        min(positive.max() + 2, len(sourcewave)))
            counts = np.trapz(sourcewave[ind] * newtrans[ind] *
                              sourceflux[..., ind],
                              sourcewave[ind], axis=-1)
            return np.squeeze(counts)
        else:
            return float('NaN')

    def obj_counts_lores(self, sourcewave, sourceflux, sourceflux_unc=0):
        """Project source spectrum onto filter and return the detector
        signal. This method differs from ``obj_counts`` in that the source
        spectrum is interpolated onto the transmission spectrum instead of
        vice-versa, which is necessary with the source spectrum does not
        adequately sample features in the transmission spectrum.

        :param sourcewave:
            Spectrum wavelength (in AA), ndarray of shape (nwave).  Must be
            monotonic increasing.

        :param sourceflux:
            Associated flux (assumed to be in erg/s/cm^2/AA), ndarray of shape
            (nspec,nwave).

        :returns counts:
            Detector signal(s) (nspec).
        """
        sourceflux = np.squeeze(sourceflux)
        assert sourceflux.ndim == 1, "Only a single source allowed."
        assert sourcewave[1] > sourcewave[0], "``sourcewave`` not in ascending order."
        # Interpolate source spectrum to filter transmission
        newflux = np.interp(self.wavelength, sourcewave, sourceflux,
                            left=0., right=0.)

        # Integrate lambda*f_lambda*R
        if True in (newflux > 0.):
            counts = np.trapz(self.wavelength * self.transmission * newflux,
                              self.wavelength)
            return np.squeeze(counts)
        else:
            return float('NaN')
    
    def ab_mag(self, sourcewave, sourceflux, lores=False, **extras):
        """Project source spectrum onto filter and return the AB magnitude.

        :param sourcewave:
            Spectrum wavelength (in AA), ndarray of shape (nwave).

        :param sourceflux:
            Associated flux (assumed to be in erg/s/cm^2/AA), ndarray of shape
            (nobj,nwave).

        :returns mag:
            AB magnitude of the source.
        """
        if not lores:
            counts = self.obj_counts(sourcewave, sourceflux, **extras)
        else:
            counts = self.obj_counts_lores(sourcewave, sourceflux, **extras)
        return -2.5 * np.log10(counts / self.ab_zero_counts)

    def vega_mag(self, sourcewave, sourceflux, lores=False, **extras):
        """Project source spectrum onto filter and return the Vega magnitude.

        :param sourcewave:
            Spectrum wavelength (in AA), ndarray of shape (nwave).

        :param sourceflux:
            Associated flux (assumed to be in erg/s/cm^2/AA), ndarray of shape
            (nobj,nwave).

        :returns mag:
            Vega magnitude of the source.
        """
        if not lores:
            counts = self.obj_counts(sourcewave, sourceflux, **extras)
        else:
            counts = self.obj_counts_lores(sourcewave, sourceflux, **extras)
        return -2.5 * np.log10(counts / self.vega_zero_counts)

# ------------
# Useful utilities
# -------------


def load_filters(filternamelist, **kwargs):
    """Given a list of filter names, this method returns a list of Filter
    objects.

    :param filternamelist:
        List of strings giving names of the filters.

    :returns filterlist:
        A list of filter objects.
    """
    return [Filter(f, **kwargs) for f in filternamelist]


def getSED(sourcewave, sourceflux, filterlist=None, **kwargs):
    """Takes wavelength vector, a flux array and list of Filter objects and
    returns the SED in AB magnitudes.

    :param sourcewave:
        Spectrum wavelength (in AA), ndarray of shape (nwave).

    :param sourceflux:
        Associated flux (assumed to be in erg/s/cm^2/AA), ndarray of shape
        (nsource,nwave).

    :param filterlist:
        List of filter objects, of length nfilt.

    :returns sed:
        array of broadband magnitudes, of shape (nsource, nfilter).
    """
    if filterlist is None:
        return None
    sourceflux = np.atleast_2d(sourceflux)
    sedshape = list(sourceflux.shape[:-1]) + [len(filterlist)]
    sed = np.zeros(sedshape)
    for i, f in enumerate(filterlist):
        sed[..., i] = f.ab_mag(sourcewave, sourceflux, **kwargs)
    return np.squeeze(sed)


def list_available_filters():
    """Return a list of filter names that are available by default, i.e. which
    have been installed in the sedpy/data/filters/ directory.
    """
    try:
        names = resource_listdir('sedpy', '/data/filters/')
    except:
        names = os.listdir(os.path.join(sedpydir, '/data/filters/'))

    parfiles = [n.replace('.par', '') for n in names if n[-4:] == '.par']
    return parfiles


def filter_dict(filterlist):
    fdict = {}
    for i, f in enumerate(filterlist):
        fdict[f.nick] = i
    return fdict

# ------------------
# Routines for spectra
# ------------------


def Lbol(wave, spec, wave_min=90, wave_max=1e6):
    """Calculate the bolometric luminosity of a spectrum or spectra.

    :param wave:
       The wavelength vector of length nwave.

    :param spec:
       The spectra, of shape (...,nsource, nwave), in F_lambda.

    :param wave_min:
       Minimum wavelength for the integral.

    :param max_wave:
       Maximum wavelength for the integral

    :returns lbol:
       The bolometric luminosity, integrated from wave_min to wave_max.  Array
       of length (...nsource)
    """
    inds = np.where(np.logical_and(wave < wave_max, wave >= wave_min))
    return np.trapz(spec[..., inds[0]], wave[inds])


def air2vac(air):
    """Convert from in-air wavelengths to vacuum wavelengths.  Based on Allen's
    Astrophysical Quantities.

    :param air:
        The in-air wavelengths.

    :returns vac:
        The corresponding vacuum wavelengths.
    """
    ss = 1e4 / air
    vac = air * (1 + 6.4328e-5 + 2.94981e-2 / (146 - ss**2) +
                 2.5540e-4 / (41 - ss**2))
    return vac


def vac2air(vac):
    """Convert from vacuum wavelengths to in-air wavelengths.  Follows the SDSS
    statement of the IAU standard from Morton 1991 ApJS.

    vac2air(air2vac(wave)) yields wave to within 1 part in a million over the
    optical range.

    :param vac:
        The vacuum wavelengths.

    :returns air:
        The corresponding in-air wavelengths.
    """
    conv = (1.0 + 2.735182e-4 + 131.4182 / vac**2 + 2.76249e8 / vac**4)
    return vac / conv


def lsf_broaden(wave, spec, lsf=None, outwave=None,
                return_kernel=False, fwhm=False, **kwargs):
    """Broaden a spectrum using a wavelength dependent line spread
    function.  This function is only approximate because it doesn't
    actually do the integration over pixels, so for sparsely sampled
    points you'll have problems.

    :param wave:
        input wavelengths
    :param lsf:
        A function that returns the gaussian dispersion at each
        wavelength.  This is assumed to be in sigma unless ``fwhm`` is
        ``True``
    :param outwave:
        Optional output wavelengths

    :param kwargs:
        Passed to lsf()

    :returns newspec:
        The broadened spectrum
    """
    if outwave is None:
        outwave = wave
    if lsf is None:
        return np.interp(outwave, wave, spec)
    dw = np.gradient(wave)
    sigma = lsf(outwave, **kwargs)
    if fwhm:
        sigma = sigma / 2.35
    kernel = outwave[:, None] - wave[None, :]
    kernel = (1 / (sigma * np.sqrt(np.pi * 2))[:, None] *
              np.exp(-kernel**2 / (2 * sigma[:, None]**2)) *
              dw[None, :])
    # should this be axis=0 or axis=1?
    kernel = kernel / kernel.sum(axis=1)[:,None]
    newspec = np.dot(kernel, spec)
    # kernel /= np.trapz(kernel, wave, axis=1)[:, None]
    # newspec = np.trapz(kernel * spec[None, :], wave, axis=1)
    if return_kernel:
        return newspec, kernel
    return newspec


def smooth_vel(wave, spec, sigma, outwave=None, inres=0):
    """
    :param sigma:
        Desired total output velocity resolution (km/s), including
        input velocity resolution
    """
    sigma_eff = np.sqrt(sigma**2 - inres**2) / 2.998e5
    if outwave is None:
        outwave = wave
    if sigma <= inres:
        if inres > 0.0:
            print("observate.smooth_vel warning: You requested a "
                  "total output velocity dispersion {0} that is lower "
                  "than the input velocity dispersion {1}.  Returning "
                  "the interpolated input spectrum".format(sigma, inres))
        return np.interp(wave, outwave, flux)

    lnwave = np.log(wave)
    flux = np.zeros(len(outwave))
    norm = 1 / np.sqrt(2 * np.pi) / sigma

    for i, w in enumerate(outwave):
        x = np.log(w) - lnwave
        f = np.exp(-0.5 * (x / sigma_eff)**2)
        flux[i] = np.trapz(f * spec, x) / np.trapz(f, x)
    return flux


def vel_broaden(sourcewave, sourceflux, sigma_in, sigma0=0,
                outwave=None, nsig=5.0, minusewave=0, maxusewave=1e8):
    """Vectorized version of velocity broadening.  This can become
    very slow when memory constraints are reached (i.e. when nw_in *
    nw_out * nsource is large).  This should be rewritten to work with
    (fft) convolutions for speed, though that requires regular (log)
    velocity scale.
    """
    sourceflux = np.atleast_2d(sourceflux)
    # Sigma after accounting for intrinsic sigma of the library.
    sigma = np.sqrt(sigma_in**2 - sigma0**2)

    if outwave is None:
        outwave = sourcewave
    # Set up arrays and limits
    minw, maxw = outwave.min(), outwave.max()
    maxw *= (1 + nsig * sigma / lightspeed)
    minw *= (1 - nsig * sigma / lightspeed)
    use = ((sourcewave > np.max([minw, minusewave])) &
           (sourcewave < np.min([maxw, maxusewave])))

    # don't need K since renormalizing weights
    # K = 1 / (sigma * np.sqrt(2 * np.pi))
    wr = outwave[:, None] / sourcewave[None, use]
    v = lightspeed / 1e13 * (1 - wr)
    ee = np.exp(-0.5 * v**2 / sigma**2)
    # renormalize - why?
    ee /= np.trapz(ee, x=v, axis=-1)[:,None]
    # now, integrate over the velocities in the sourcewave direction
    flux = np.trapz(sourceflux[:, None, use] * ee[None, :, :],
                    x=v[None, :, :], axis=-1)

    # H = lightspeed/1e13/sigma
    # flux = np.trapz( sourceflux[:,None,use] *
    #                  np.exp(-0.5 * (H * (1 - outwave[:,None] /
    #                                      sourcewave[None, use ]))**2),
    #                  x= lightspeed/1e13 * (1 - outwave[None,:,None]/
    #                                        sourcewave[None,None,:]), axis = -1)

    return flux  # * K


def vel_broaden_fast(sourcewave, sourceflux, sigma):
    pass


def smooth_wave(wave, spec, sigma, outwave=None,
                inres=0, in_vel=False, **extras):
    """
    :param sigma:
        Desired reolution in wavelength units

    :param inres:
        Resolution of the input, in either wavelength units or
        lambda/dlambda (c/v)

    :param in_vel:
        If True, the input spectrum has been smoothed in velocity
        space, and inres is in dlambda/lambda.
    """
    if outwave is None:
        outwave = wave
    if inres <= 0:
        sigma_eff = sigma
    elif in_vel:
        sigma_min = np.max(outwave) / inres
        if sigma < sigma_min:
            raise ValueError("Desired wavelength sigma is lower "
                             "than the value possible for this input "
                             "spectrum ({0}).".format(sigma_min))
        # Make an approximate correction for the intrinsic wavelength
        # dependent dispersion.  This doesn't really work.
        sigma_eff = np.sqrt(sigma**2 - (wave / inres)**2)
    else:
        if sigma < inres:
            raise ValueError("Desired wavelength sigma is lower "
                             "than the value possible for this input "
                             "spectrum ({0}).".format(sigma_min))
        sigma_eff = np.sqrt(sigma**2 - inres**2)

    flux = np.zeros(len(outwave))
    for i, w in enumerate(outwave):
        # if in_vel:
        #     sigma_eff = np.sqrt(sigma**2 - (w/inres)**2)
        x = (wave - w) / sigma_eff
        f = np.exp(-0.5 * x**2)
        flux[i] = np.trapz(f * spec, wave) / np.trapz(f, wave)
    return flux


def wave_broaden(sourcewave, sourceflux, fwhm, fwhm0=0, outwave=None,
                 nsig=5.0, minusewave=0, maxusewave=1e8):
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
    maxw *= (1 + nsig * sigma)
    minw *= (1 - nsig * sigma)
    use = ((sourcewave > np.max([minw, minusewave])) &
           (sourcewave < np.min([maxw, maxusewave])))

    dl = outwave[:,None] - sourcewave[None, use]
    ee = np.exp(-0.5 * dl**2 / sigma**2)
    ee /= np.trapz(ee, x=dl, axis=-1)[:,None]
    flux = np.trapz(sourceflux[:, None, use] * ee[None, :, :],
                    x=dl[None, :, :], axis=-1)

    return flux


def broaden(sourcewave, sourceflux, width, width0=0,
            stype='vel', **kwargs):
    if stype == 'vel':
        return vel_broaden(sourcewave, sourceflux, width,
                           sigma0=width0, **kwargs)
    elif stype == 'wave':
        return wave_broaden(sourcewave, sourceflux, width,
                            fwhm0=width0, **kwargs)

#    K=(sigma*sqrt(2.*!PI))^(-1.)
#    for iw in range(len(outwave)):
#        dl=(outwave[iw]-wave)**2.
#        this=np.where(dl < length*sigma**2.,count)
#        if count > 0:
#            ee=exp(-0.5*dl[this]/sigma^2)
#            broadflux[iw]=int_tabulated(wave[this],flux[this]*ee)
#    return,broadflux

# def redshift(sourcewave,sourceflux, z):


def smooth(wave, spec, sigma, smooth_velocity=True, **kwargs):
    if smooth_velocity:
        return smooth_vel(wave, spec, sigma, **kwargs)
    else:
        return smooth_wave(wave, spec, sigma, **kwargs)


def selftest():
    """Compare to the values obtained from the K-correct code
    (which uses a slightly different Vega spectrum)
    """
    filternames = ['galex_FUV',
                   'sdss_u0','sdss_g0','sdss_r0','sdss_i0',
                   'spitzer_irac_ch2']
    weff_kcorr = [1528.0,
                  3546.0, 4669.6, 6156.25, 7471.57,
                  44826.]
    msun_kcorr = [18.8462,
                  6.38989, 5.12388, 4.64505, 4.53257,
                  6.56205]
    ab2vega_kcorr = [2.3457,
                     0.932765, -0.0857, 0.155485, 0.369598,
                     3.2687]

    filterlist = loadFilters(filternames)
    for i in range(len(filterlist)):
        print(filterlist[i].wave_effective, filterlist[i].solar_ab_mag,
              filterlist[i].ab_to_vega)
        assert (abs(filterlist[i].wave_effective - weff_kcorr[i]) <
                (weff_kcorr[i] * 0.01))
        assert abs(filterlist[i].solar_ab_mag - msun_kcorr[i]) < 0.05
        # the below fails because of the vega spectrum used by k_correct
        # assert abs(filterlist[i].ab_to_vega+ab2vega_kcorr[i]) < 0.05


def selftest_broaden():
    import fsps
    sps = fsps.StellarPopulation()
    ages = [1., 2., 3., 4., 5., 6., 7.]
    allspec = np.zeros([len(ages), len(sps.wavelengths)])
    for i,tage in enumerate(ages):
        wave, spec = sps.get_spectrum(peraa=True, tage=tage)
        allspec[i,:] = spec

    outwave = wave[(wave > 1e3) & (wave < 1e4)]
    vbflux = vel_broaden(wave, allspec, 150., outwave=outwave)
    wbflux = wave_broaden(wave, allspec, 5., outwave=outwave)
