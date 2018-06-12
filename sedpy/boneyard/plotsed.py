# This code is an unfinished attempt to plot SEDs as probability
# distributions p(flux, wave)

import numpy as np

class Bunch(object):
    """ Simple storage.
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class SedPoint(object):

    def __init__(self, flux=0., unc=0., band=None):

        self.flux = flux
        self.unc = unc
        self.band = band

    @property
    def var(self):
        return self.unc**2

    @property
    def pnorm(self):
        return (self.var * 2 * np.pi)**(-0.5)

    @property
    def wave_min(self):
        return self.band.wave_effective - self.band.effective_width/2.

    @property
    def wave_max(self):
        return self.band.wave_effective + self.band.effective_width/2.

    def flux_lnprob(self, fluxes):
        fnp =  (np.array(fluxes) - self.flux)**2 / (2 * self.var) 
        return  np.log(self.pnorm) - fnp

    def spans(self, wavelengths):
        waves = np.array(wavelengths)
        return (waves < self.wave_max) & (waves > self.wave_min)

    def sed_prob(self, wavelength_grid, flux_grid):
        return self.spans(wavelength_grid), self.flux_lnprob(flux_grid)


def sed_to_psed(filters, fluxes, uncs, wgrid, fgrid):
    sed = []
    lnpsed = np.zeros(( len(wgrid), len(fgrid)))
    for filt, flux, unc in zip(filters, fluxes, uncs):
        sed_point =  SedPoint(flux, unc, filt)
        winds, lnprob = sed_point.sed_prob(wgrid, fgrid)
        lnpsed[winds,:] += lnprob[None, :]
        sed += [sed_point]
        lnpsed -= np.log(np.trapz(np.exp(lnpsed), fgrid, axis = 1))[:, None]
        lnpsed += np.log(np.trapz(np.ones(len(fgrid)), fgrid))
    return lnpsed, sed


def test():
    from sedpy import observate
    import fsps
    import matplotlib.pyplot as pl

    filters = ['galex_NUV', 'sdss_u0', 'sdss_r0', 'sdss_r0', 'sdss_i0', 'sdss_z0',
               'bessell_U', 'bessell_B', 'bessell_V', 'bessell_R', 'bessell_I',
               'twomass_J','twomass_H']
    flist = observate.load_filters(filters)

    sps = fsps.StellarPopulation(compute_vega_mags=False)
    wave, spec = sps.get_spectrum(tage=1.0, zmet=2, peraa=True)

    sed = observate.getSED(wave, spec, flist)
    sed_unc = np.abs(np.random.normal(1, 0.3, len(sed)))

    wgrid = np.linspace( 2e3, 13e3, 1000)
    fgrid = np.linspace( -13, -9, 100)
    psed, sedpoints = sed_to_psed(flist, sed, sed_unc, wgrid, fgrid)

    pl.imshow(np.exp(psed).T, cmap='Greys_r',
              interpolation='nearest', origin ='upper', aspect='auto')
    
