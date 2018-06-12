import numpy as np


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
