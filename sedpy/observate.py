# -*- coding: utf-8 -*-

"""observate.py - Python module for storing filter information and tools for
projecting spectra onto filters.  Also includes tools for convolving spectra.

Assumed input units are erg/s/cm^2/AA and AA
"""

import numpy as np
import os
try:
    from pkg_resources import resource_filename, resource_listdir
except(ImportError):
    pass
from .reference_spectra import vega, solar, sedpydir


__all__ = ["Filter", "FilterSet",
           "load_filters", "list_available_filters", "getSED",
           "air2vac", "vac2air", "Lbol"]


lightspeed = 2.998e18  # AA/s


class Filter(object):
    """This class operates on filter transmission files.  It reads simple
    2-column ascii files containing filter transmission curves. It caches a
    number of useful filter quantities.  Methods are provided to project a
    source spectrum onto the filter and return the magnitude, AB or Vega.

    Parameters
    ----------
    kname : string (default: 'sdss_r0')
        The kcorrect style name of the filter, excluing '.par',
        e.g. sdss_r0.

    nick : string (optional)
        A nickname to associate with the filter.

    directory : string (optional)
        The path to the directory containing the filter file.  If not given
        then the sedpy/data/filters/ directory will be searched.

    trans_colnum : int or None (default, None)
        Column from which to pull the transmission, starting with 1. Will attempt
        to update the filter name with the column name.

    data: None or tuple of ndarray (optional, default None)
        If provided, a 2-tuple of ndarrays giving wavelength (AA) and
        transmission (detector signal per photon).

    dlnlam : float (optional)
        If given, interpolate the transmission curve onto a grid in
        :math:`\ln(\lambda)` with this spacing (and begining at ``wmin``)

    wmin : float (optional, default 100.)
        If ``dlnlam`` is supplied then ``wmin`` gives the minimum wavelength
        (in Angstroms) for the grid.

    min_trans : float (optional, default 1e-5)
        The minimum transmission value (as fraction of the maximum transmission
        value) to consider as a useful positive value.  Useful for curves that
        have extremely low amplitude components far from the main part of the
        filter when speed is more important than accuracy.
    """
    ab_gnu = 3.631e-20  # AB reference spctrum in erg/s/cm^2/Hz
    npts = 0

    def __init__(self, kname='sdss_r0', nick=None,
                 directory=None, trans_colname=None, data=None,
                 dlnlam=None, wmin=1e2, min_trans=1e-5,  **extras):
        """Read in the filter data, calculate and cache some properties of the
        filters.
        """
        self.name = kname
        if nick is None:
            self.nick = kname
        else:
            self.nick = nick

        self.min_trans = min_trans
        self.trans_colname = trans_colname

        if directory is None:
            loc = os.path.join('data', 'filters', kname + '.par')
            try:
                self.filename = resource_filename("sedpy", loc)
            except:
                self.filename = os.path.join(sedpydir, loc)
        else:
            self.filename = os.path.join(directory, kname+'.par')

        # Load wave, trans directly from data arrays
        if data is not None:
            wave, trans = data
            self._process_filter_data(wave, trans)
            self.filename = None

        if isinstance(self.filename, str):
            if not os.path.isfile(self.filename):
                raise ValueError('Filter transmission file {0} does '
                                 'not exist!'.format(self.filename))
            self.load_filter(self.filename)
        elif isinstance(self.filename, type(None)):
            pass
        else:
            msg = "{} of type {} is not a valid transmission file name: check your string types?"
            raise TypeError(msg.format(self.filename, type(self.filename)))

        if dlnlam is not None:
            self.gridify_transmission(dlnlam, wmin=wmin)

        self.get_properties()

    def __repr__(self):
        return '{}({})'.format(self.__class__, self.name)

    def load_filter(self, filename):
        """Read a filter in simple two column ascii format and populate the
        wavelength and transmission arrays.  The first column is wavelength in
        AA and the second column is transmission (detector signal per photon)

        Parameters
        ----------
        filename : string
            The fully qualified path and filename of the file that contains the
            filter transmission.
        """
        if self.trans_colname is None:
            # basic
            wave, trans = np.genfromtxt(filename, usecols=(0, 1), unpack=True)
        else:
            # get colnum from colname
            with open(filename, "r") as f:
                line = "#"
                while (line[0] in "#\n"):
                    cols = line[1:].split()
                    line = f.readline()
            try:
                colnum = cols.index(self.trans_colname)
            except:
                print(f"{self.trans_colname} not in {cols}")
                raise
            wave, trans = np.genfromtxt(filename, usecols=(0, colnum), unpack=True)

            self.name = f"{self.name}_{self.trans_colname}"
            self.nick = f"{self.nick}_{self.trans_colname}"

        # Clean negatives, NaNs, and Infs, then sort, then store
        self._process_filter_data(wave, trans)

    @property
    def transmission(self):
        return self._transmission

    @property
    def wavelength(self):
        return self._wavelength

    def _process_filter_data(self, wave, trans):
        """Clean up transmission data and assign to attributes.

        Parameters
        ----------
        wave : ndarray of shape ``(N_wave)``
            Wavelength, in Angstroms, of the transmission curve

        trans : ndarray of shape ``(N_wave)``
            Filter transmission at the wavelengths corresponding to ``wave``
        """
        ind = np.isfinite(trans) & (trans >= 0.0)
        order = wave[ind].argsort()
        self.npts = ind.sum()
        self._wavelength = wave[ind][order]
        self._transmission = trans[ind][order]
        self._remove_extra_zeros(self.min_trans)

    def _remove_extra_zeros(self, min_trans=0):
        """Remove extra leading or trailing zero transmission points.  This
        leaves one zero before (after) the first (last) non-zero transmission
        point, if they were present in the original transmission function.

        Parameters
        ----------
        min_trans : float (optional, default: 1e-5)
            Defines zero, in terms of fraction of the maximum transmission.
        """
        v = np.argwhere(self.transmission > (self.transmission.max() * min_trans))
        inds = slice(max(v.min() - 1, 0), min(v.max() + 2, len(self.transmission)))
        self._wavelength = self._wavelength[inds]
        self._transmission = self._transmission[inds]
        self.npts = len(self.wavelength)

    def gridify_transmission(self, dlnlam, wmin=1e2):
        """Place the transmission function on a regular grid in lnlam
        (angstroms) defined by a lam_min and dlnlam.  Note that only the
        non-zero values of the transmission on this grid stored.  The indices
        corresponding to these values are stored as the `inds` attribute (a
        slice object). (with possibly a zero at either end.)

        Parameters
        ----------
        dlnlam : float
            The spacing in ln-lambda of the regular wavelength grid onto which
            the filter is to be placed.

        wmin : float (optional, default: 100)
            The starting wavelength (Angstroms) for the regular grid.
        """
        # find min and max of the filter in a global wavelength grid given by
        # wmin, wmax, and dlnlam
        ind_min = int(np.floor((np.log(self.wavelength.min()) - np.log(wmin)) / dlnlam))
        ind_max = int(np.ceil((np.log(self.wavelength.max()) - np.log(wmin)) / dlnlam))
        lnlam = np.linspace(ind_min * dlnlam + np.log(wmin),
                            ind_max * dlnlam + np.log(wmin), ind_max - ind_min)
        lam = np.exp(lnlam)
        # TODO: replace with a rebinning
        trans = np.interp(lam, self.wavelength, self.transmission,
                          left=0., right=0.)

        self.wmin = wmin
        self.dlnlam = dlnlam
        self.inds = slice(ind_min, ind_max)
        self._wavelength = lam
        self._transmission = trans
        self.dwave = np.gradient(lam)

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

        tmax = self.transmission.max()
        self.wave_effective = np.exp(i0 / i1)
        self.wave_pivot = np.sqrt(i2 / i1)
        self.wave_mean = self.wave_effective
        self.wave_average = i2 / i3
        self.rectangular_width = i3 / tmax

        wpeak = self.wavelength[np.argmax(self.transmission)]
        sel = (self.wavelength < wpeak) & (self.transmission < 0.5 * tmax)
        bind = np.argmax(self.wavelength[sel])
        self.blue_edge = self.wavelength[sel][bind]
        sel = (self.wavelength > wpeak) & (self.transmission < 0.5 * tmax)
        rind = np.argmin(self.wavelength[sel])
        self.red_edge = self.wavelength[sel][rind]

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
            self.vega_zero_counts = self.obj_counts(vega[:, 0], vega[:, 1])
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
        """Plot the filter transmission curve
        """
        if self.npts > 0:
            if ax is None:
                import matplotlib.pyplot as pl
                _, ax = pl.subplots()
                ax.set_title(self.nick)
            if normalize:
                ax.plot(self.wavelength, self.transmission / self.transmission.max())
            else:
                ax.plot(self.wavelength, self.transmission)
            return ax

    def obj_counts_hires(self, sourcewave, sourceflux, sourceflux_unc=0):
        """Project source spectrum onto filter and return the detector signal.

        Parameters
        ----------
        sourcewave : ndarray of shape ``(N_pix,)``
            Spectrum wavelength (in Angstroms). Must be monotonic increasing.

        sourceflux : ndarray of shape ``(N_source, N_pix)``
            Associated flux (assumed to be in erg/s/cm^2/AA)

        Returns
        -------
        counts : ndarray of shape ``(N_source,)``
            Detector signal(s).
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
        signal. This method differs from ``obj_counts_hires`` in that the
        source spectrum is interpolated onto the transmission spectrum instead
        of vice-versa, which is necessary with the source spectrum does not
        adequately sample features in the transmission spectrum.

        Parameters
        ----------
        sourcewave : ndarray of shape ``(N_pix,)``
            Spectrum wavelength (in Angstroms). Must be monotonic increasing.

        sourceflux : ndarray of shape ``(N_source, N_pix)``
            Associated flux (assumed to be in erg/s/cm^2/AA)

        Returns
        -------
        counts : ndarray of shape ``(N_source,)``
            Detector signal(s).
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

    def obj_counts_grid(self, sourceflux, source_offset=None):
        """Project source spectrum onto filter and return the detector
        signal. This method differs from ``obj_counts_*res`` in that the source
        spectrum is assumed to be interpolated onto a logarithmic grid in
        lambda with spacing and minimum wavelength given by `dlnlam` and `wmin`
        filter attributes.

        Parameters
        ----------
        sourceflux : ndarray of shape ``(N_source, N_pix)``
            Source flux (assumed to be in erg/s/cm^2/AA) assumed to be on the
            logarithmic grid defined by ``wmin`` and ``dlnlam``, with wavelength
            ascending.

        source_offset : int (optional)
            The (hypothetical) element of the sourceflux array corresponding to
            observed frame wavelength of ``wmin``.  Should be negative if the first
            element of sourcewave corresponds to a wavelength > wmin.

        Returns
        -------
        counts : ndarray of shape ``(N_source,)``
            Detector signal(s).
        """
        # Note that redshifting can also be incorporated using negative `source_offset`,
        # modulo factors of (1+z).  I think.

        if source_offset is None:
            valid = self.inds
        else:
            valid = slice(self.inds.start + source_offset,
                          self.inds.stop + source_offset)
        counts = np.sum(sourceflux[..., valid] * self.transmission *
                        self.wavelength * self.dwave, axis=-1)
        return np.squeeze(counts)

    def obj_counts(self, sourcewave, sourceflux, lores=False, gridded=False, **extras):
        """Project a spectrum onto a filter and return the detector signal.
        This method uses the keywords `lores` and `gridded` to choose between
        the various projection algorithms.

        Parameters
        ----------
        sourcewave : ndarray of shape ``(N_pix,)``
            Spectrum wavelength (in Angstroms). Must be monotonic increasing.

        sourceflux : ndarray of shape ``(N_source, N_pix)``
            Associated flux (assumed to be in erg/s/cm^2/AA)

        lores : bool (optional, default: False)
            Switch to interpolate the source spectrum onto the wavelength grid
            of the filter transmission curve, instead of the vice-versa (which
            is the default).  Useful for a spectrum with sparse sampling in
            wavelength compared to the filter transmission curve (e.g. narrow
            bands and BaSeL lores spectra)

        gridded : bool (optional, default: False)
            Switch to accomplish the filter projection via simple sums.  This
            can be faster, but assumes the sourcewave and sourceflux are on a
            logarithmic wavelength grid given by the dlnlam and wmin attributes
            of the Filter object.

        Returns
        -------
        counts : ndarray of shape ``(N_source,)``
            Detector signal(s).

        """
        if gridded:
            counts = self.obj_counts_grid(sourceflux, **extras)
        elif lores:
            counts = self.obj_counts_lores(sourcewave, sourceflux, **extras)
        else:
            counts = self.obj_counts_hires(sourcewave, sourceflux, **extras)

        return counts

    def ab_mag(self, sourcewave, sourceflux, **extras):
        """Project source spectrum onto filter and return the AB magnitude.

        Parameters
        ----------
        sourcewave : ndarray of shape ``(N_pix,)``
            Spectrum wavelength (in Angstroms). Must be monotonic increasing.

        sourceflux : ndarray of shape ``(N_source, N_pix)``
            Associated flux (assumed to be in erg/s/cm^2/AA)

        Returns
        -------
        mag : float or ndarray of shape ``(N_source,)``
            AB magnitude of the source(s).
        """
        counts = self.obj_counts(sourcewave, sourceflux, **extras)
        return -2.5 * np.log10(counts / self.ab_zero_counts)

    def vega_mag(self, sourcewave, sourceflux, **extras):
        """Project source spectrum onto filter and return the Vega magnitude.

        Parameters
        ----------
        sourcewave : ndarray of shape ``(N_pix,)``
            Spectrum wavelength (in Angstroms). Must be monotonic increasing.

        sourceflux : ndarray of shape ``(N_source, N_pix)``
            Associated flux (assumed to be in erg/s/cm^2/AA)

        Returns
        -------
        mag : float or ndarray of shape ``(N_source,)``
            Vega magnitude of the source(s).
        """
        counts = self.obj_counts(sourcewave, sourceflux, **extras)
        return -2.5 * np.log10(counts / self.vega_zero_counts)


class FilterSet(object):
    """
    A class representing a set of filters, and enabling fast projections of
    sources onto these filters via pre-computation, caching, and a global
    wavelength array.
    """

    def __init__(self, filterlist, wmin=None, dlnlam=None,
                 **loading_kwargs):
        """
        Parameters
        ----------
        filterlist : list of strings
            The filters to include in this set, in order

        wmin : float
            Minimum wavelength (AA) for the wavelength grid for the FilterSet

        dlnlam : float
            logarithmic spacing for the wavelength grid of the FilterSet
        """
        self.filternames = filterlist
        native = load_filters(self.filternames, **loading_kwargs)
        self._set_filters(native, wmin=wmin, dlnlam=dlnlam, **loading_kwargs)
        self._build_super_trans()

    def __len__(self):
        return len(self.filternames)

    def _set_filters(self, native, wmin=None, wmax=None, dlnlam=None,
                     **loading_kwargs):
        """Set the filters and the wavelength grid

        native : list of Filter() instances
            The Filter objects that will be part of this FilterSet, including
            their native transmission
        """
        # get the wavelength grid parameters
        self.dlnlam_native = np.array([np.diff(np.log(f.wavelength)).min()
                                       for f in native])
        if dlnlam is None:
            dlnlam = min(self.dlnlam_native.min(), 1e-3)
        if wmin is None:
            wmin = np.min([f.wavelength.min() for f in native])
        if wmax is None:
            wmax = np.max([f.wavelength.max() for f in native])
        self.wmin = wmin
        self.wmax = np.exp(np.log(wmax) + dlnlam)
        self.dlnlam = dlnlam

        # Build the wavelength grid
        self.lnlam = np.arange(np.log(self.wmin), np.log(self.wmax), self.dlnlam)
        self.lam = np.exp(self.lnlam)

        # build the filters on that grid
        self.filters = load_filters(self.filternames, wmin=self.wmin,
                                    dlnlam=self.dlnlam, **loading_kwargs)

    def _build_super_trans(self):
        """Build the global gridded transmission array. Populates FilterSet.trans, which
        is an array of shape (nfilters, nlam) giving :math:`T\times \lambda \times
        d\lambda / Z` where `T` is the detector signal per photon and `Z` is the
        detector signal for a 1 maggie source
        """
        # TODO: investigate using arrays that cover different wavelength
        # intervals (but same number of elements) for each filter. Then when
        # integrating (with dot or vectorized trapz) the source spectrum must be
        # offset for each filter. This would remove the vast majority of zeros

        ab_counts = np.array([f.ab_zero_counts for f in self.filters])

        self.trans = np.zeros([len(self.filters), len(self.lam)])
        for j, f in enumerate(self.filters):
            dwave = np.gradient(self.lam[f.inds])
            self.trans[j, f.inds] = f.transmission[:] * f.wavelength[:] * dwave / ab_counts[j]

        # this is for the numba unequal inmplementation
        self.trans_list = [self.trans[j, f.inds] for j, f in enumerate(self.filters)]
        self.frange = np.array([(f.inds.start, f.inds.stop) for f in self.filters])

    def interp_source(self, inwave, sourceflux):
        """Interpolate spectra onto the global wavelength grid.

        Parameters
        ----------
        inwave : ndarray of shape ``(N_pix,)``
            Input wavelength (Angstroms)

        sourceflux : ndarray of shape ``(N_source, N_pix)``
            Input source fluxes at the wavelengths given by `inwave`

        Returns
        -------
        flux : ndarray of shape ``(N_pix,)`` or ``(N_pix, N_source)``
           Interpolated flux. Note the transpose relative to input if ``N_source > 1``
        """
        sourceflux = np.atleast_2d(sourceflux)
        # TODO: replace with scipy interpolate for vectorization?
        # TODO replace with rebinning!!!!
        flux = np.array([np.interp(self.lam, inwave, s, left=0, right=0)
                         for s in sourceflux])
        return np.squeeze(flux.T)

    def get_sed_maggies(self, sourceflux, sourcewave=None):
        """Compute the flux in maggies of the source through the given filters.

        Parameters
        ----------
        sourceflux : ndarray of shape ``(N_source, N_pix)``
            Source flux (assumed to be in erg/s/cm^2/AA)

        sourcewave: ndarray (optional)
            If supplied, the wavelength in Angstroms.  Otherwise it is assumed
            to be on the same wavelength grid as FilterSet.lam, saving the cost
            of an interpolation.

        Returns
        -------
        maggies : ndarray of shape ``(N_filter,)`` or ``(N_source, N_filter)``
            The flux (in maggies) of the sources through the filters.
        """
        if sourcewave is not None:
            flux = self.interp_source(sourcewave, sourceflux)
        else:
            flux = sourceflux.T
            assert len(flux) == len(self.lam)
        maggies = np.dot(self.trans, flux)

        return maggies.T

    def _build_trans_chunks(self):
        # NOTE: For USE WITH vectorized TRAPZ
        ab_counts = np.array([f.ab_zero_counts for f in self.filters])
        self.fstart = np.array([f.inds.start for f in self.filters])
        self.fstop = np.array([f.inds.stop for f in self.filters])
        self.nf = self.fstop - self.fstart

        # build an index array in lnlam that produces (nfilter, nf.max()) shape array
        self.finds = self.fstart[:, None] + range(self.nf.max())
        finds_right = self.fstop[:, None] - range(self.nf.max())[::-1]
        right = self.finds.max(axis=-1) >= len(self.lam)
        self.finds[right, :] = finds_right[right, :]

        # fill the transmission arrays, shifting some to the right (i.e. zero padding on the left)
        # to make them fit.  Note this is T \times \lambda, but *not* times dw
        self.trans_dw = np.zeros([len(self.filters), self.nf.max()])
        for j, f in enumerate(self.filters):
            #dwave = np.gradient(self.lam[f.inds])
            t = f.transmission[:] * f.wavelength[:] / ab_counts[j]
            if right[j]:
                self.trans_dw[j, -self.nf[j]:] = t
            else:
                self.trans_dw[j, :self.nf[j]] = t

    def get_sed_maggies_trapz(self, sourceflux, sourcewave=None):
        """Use vectorized trapz to project filters.  This can be
        slower than the default get_sed_maggies
        """
        # NOTE: THIS IS SLOWER THAN GET_SED_MAGGIES
        if sourcewave is not None:
            flux = self.interp_source(sourcewave, sourceflux)
        else:
            assert len(sourceflux) == len(self.lam)
            flux = sourceflux
        y = self.trans_dw * flux[self.finds]
        x = self.lam[self.finds]
        maggies = np.trapz(y, x, axis=-1)
        return maggies

    def get_sed_maggies_numba(self, sourceflux, sourcewave=None):
        if sourcewave is not None:
            source = self.interp_source(sourcewave, sourceflux)
        else:
            source = sourceflux

        # Allow for uneven transmission curve vectors
        maggies = _project_filters_fast(len(self), self.frange, self.trans_list, source)

        return maggies

try:
    import numba
    @numba.njit
    def _project_filters_fast(nfilt, filter_inds,
                              trans_list, source_vector):
        maggies = np.empty(nfilt)
        for i in range(nfilt):
            start, stop = filter_inds[i]
            source = source_vector[start:stop]
            maggies[i] = np.sum(trans_list[i] * source)
        return maggies

    #assert _project_filters_fast.nopython_signatures

except(ImportError):
    def _project_filters_fast(nfilt, filter_inds,
                              trans_list, source_vector):
        raise ImportError("Could not import numba")

# ------------
# Useful utilities
# -------------


def rebin(outwave, wave, trans):
    """Rebin (instead of interpolate) transmission onto given wavelength array,
    preserving total transmission.  This is useful if the output wavelength
    array is more coarseley sampled than the native transmission array, to avoid
    interpolating over important sharp features in the transmission function

    Stripped down version of specutils.FluxConservingResampler

    FIXME: implement as NUFFT?  Do something about the jaggedness (high frequency noise)
    introduced when the output is more finely sampled than the input?
    """
    assert np.all(np.diff(wave) > 0)
    assert np.all(np.diff(outwave) > 0)

    # Assume input is based on centers, so we need to work out edges
    edges = (wave[:-1] + wave[1:]) / 2
    edges = np.concatenate([[2*edges[0]-edges[1]], edges, [2*edges[-1] - edges[-2]]])
    inlo, inhi = edges[:-1], edges[1:]
    edges = (outwave[:-1] + outwave[1:]) / 2
    edges = np.concatenate([[2*edges[0]-edges[1]], edges, [2*edges[-1] - edges[-2]]])
    outlo, outhi = edges[:-1], edges[1:]

    # Make a huge matrix of the contribution of each input bin to each output bin
    # Basically a brute force convolution with variable width square kernel
    l_inf = np.where(inlo > outlo[:, None], inlo, outlo[:, None])
    l_sup = np.where(inhi < outhi[:, None], inhi, outhi[:, None])
    resamp_mat = (l_sup - l_inf).clip(0) * (inhi - inlo)

    # remove any contribution to output bins that aren't totally within the interval of the original bins.
    left_clip = np.where(outlo - inlo[0] < 0, 0, 1)
    right_clip = np.where(inhi[-1] - outhi < 0, 0, 1)
    keep_overlapping_matrix = left_clip * right_clip
    resamp_mat *= keep_overlapping_matrix[:, None]

    # Do the convolution
    out = np.dot(resamp_mat, trans) / resamp_mat.sum(axis=-1)

    return out


def load_filters(filternamelist, **kwargs):
    """Given a list of filter names, this method returns a list of Filter
    objects.

    Parameters
    ----------
    filternamelist : list of strings
        The names of the filters.

    Returns
    -------
    filterlist : list of Filter() instances
        A list of filter objects.
    """
    return [Filter(f, **kwargs) for f in filternamelist]


def getSED(sourcewave, sourceflux, filterlist=None, linear_flux=False, **kwargs):
    """Takes wavelength vector, a flux array and list of Filter objects and
    returns the SED in AB magnitudes.

    Parameters
    ----------
    sourcewave : ndarray of shape ``(N_pix,)``.
        Spectrum wavelength (in Angstroms)

    sourceflux : ndarray of shape ``(N_pix,)`` or ``(N_source, N_pix)``
        Associated flux (assumed to be in erg/s/cm^2/AA)

    filterlist : list of Filter() instances or a FilterSet() instance.
        The filters onto which to project the spectra

    linear_flux : bool, optional (default: True)
        If True return maggies, else AB mags

    Returns
    -------
    sed : ndarray of shape ``(N_source, N_filter)`` or ``(N_filter,)``
        broadband magnitudes or maggies.
    """
    if hasattr(filterlist, "get_sed_maggies"):
        maggies = filterlist.get_sed_maggies(sourceflux, sourcewave=sourcewave)
        if linear_flux:
            return maggies
        else:
            return -2.5 * np.log10(maggies)

    elif filterlist is None:
        return None

    # standard loop over filters
    sflux = np.squeeze(sourceflux)
    sedshape = sflux.shape[:-1] + (len(filterlist),)
    sed = np.zeros(sedshape)
    for i, f in enumerate(filterlist):
        sed[..., i] = f.ab_mag(sourcewave, sflux, **kwargs)
    #sed = np.atleast_1d(np.squeeze(sed))
    if linear_flux:
        return 10**(-0.4 * sed)
    else:
        return sed


def list_available_filters():
    """Return a list of filter names that are available by default, i.e. which
    have been installed in the sedpy/data/filters/ directory.
    """
    try:
        names = resource_listdir('sedpy', '/data/filters/')
    except:
        names = os.listdir(os.path.join(sedpydir, '/data/filters/'))

    parfiles = [n.replace('.par', '') for n in names if n[-4:] == '.par']
    parfiles.sort()
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

    Parameters
    ----------
    wave : ndarray
       The wavelength vector of length nwave.

    spec : ndarray
       The spectra, of shape (...,nsource, nwave), in F_lambda.

    wave_min : float
       Minimum wavelength for the integral.

    max_wave : float
       Maximum wavelength for the integral

    Returns
    -------
    lbol :
       The bolometric luminosity, integrated from wave_min to wave_max.  Array
       of length (...nsource)
    """
    inds = np.where(np.logical_and(wave < wave_max, wave >= wave_min))
    return np.trapz(spec[..., inds[0]], wave[inds])


def air2vac(air):
    """Convert from in-air wavelengths to vacuum wavelengths.  Based on Allen's
    Astrophysical Quantities.

    Parameters
    ----------
    air : ndarray of shape ``(N_pix,)``
        The in-air wavelengths.

    Returns
    -------
    vac: ndarray of shape ``(N_pix,)``
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

    Parameters
    ----------
    vac: ndarray of shape ``(N_pix,)``
        The vacuum wavelengths.

    Returns
    -------
    air : ndarray of shape ``(N_pix,)``
        The corresponding in-air wavelengths.
    """
    conv = (1.0 + 2.735182e-4 + 131.4182 / vac**2 + 2.76249e8 / vac**4)
    return vac / conv
