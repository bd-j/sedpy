# -*- coding: utf-8 -*-

"""modelgrid.py  - Deprecated code for interpolating regular or irregular grids
of models.  Use `scipy.interpolate.RegularGridInterpolator` or
`scipy.interpolate.LinearNDInterpolator`.

Kept here since it is used by some other (old) packages of mine.
"""

import numpy as np
from scipy.spatial import Delaunay
try:
    import astropy.io.fits as pyfits
except (ImportError):
    import pyfits
try:
    import sklearn.neighbors
except(ImportError):
    pass
try:
    from . import observate
except (ImportError):
    print('Warning - observate not imported, SpecLibrary class unavailable')


class ModelLibrary(object):
    """Class to deal with (irregular) grids of models.  Primary
    attribute is `pars`: a structured array of parameter values of
    shape (ngrid).  Methods are provided to manipulate the parameter
    structure, and to obtain interpolation weights for specific points
    within the grid using a variety of methods.
    """
    # force generation of the DT and KDtree in first pass
    triangle_dirtiness = 1

    def __init__(self, pars=None, parnames=None):
        if pars is not None:
            self.set_pars(pars, parnames)

    def set_pars(self, pars, parnames):
        self.pars = self.structure_array(pars, parnames)
        self.ngrid = self.pars.shape[0]

    def add_par(self, value, name, dtype='<f8'):
        newp = self.structure_array(value, [name], types=[dtype])
        self.pars = self.join_struct_arrays([self.pars, newp])
        pass

    def par_names(self):
        return self.pars.dtype.names

    def par_range(self, parname, inds=None):
        range_list = [(np.nanmin(self.pars[inds][p]),
                       np.nanmax(self.pars[inds][p]), p)
                      for p in parname]
        return range_list

    def structure_array(self, values, fieldnames, types=None):
        """Turn a numpy array of floats into a structurd array.

        :param values: ndarray, shape(nrec, nfield).
            Values of the structure.

        :param fieldnames: string sequence of length nfield.
            A list or string array of parameter names with length
            NFIELD.

        :param types: string sequence, optional
            Format identifiers for each field. Defaults to '<f8'.

        :returns struct:
            A numpy structured array.
        """
        # shape coercion
        values = np.atleast_2d(values)
        if values.shape[-1] != len(fieldnames):
            if values.shape[0] == len(fieldnames):
                values = values.T
            else:
                raise ValueError('array and fieldnames do not have consistent shapes!')
        nobj = values.shape[0]

        # create the dtype and structured array
        if types is None:
            types = ['<f8'] * len(fieldnames)
        dt = np.dtype(zip(fieldnames, types))
        struct = np.zeros(nobj, dtype=dt)
        for i, f in enumerate(fieldnames):
            struct[f] = values[..., i]
        return struct

    def join_struct_arrays(self, arrays):
        """From some dudes on StackOverflow.  Add equal length
        structured arrays to produce a single structure with fields
        from both.

        :param arrays:
            A sequence of structured arrays.  These must each have the
            same length.

        :returns newarray:
            A single array containing the fields of ``arrays``.
        """
        if False in [len(a) == len(arrays[0]) for a in arrays]:
            raise ValueError('array lengths do not match.')

        newdtype = np.dtype(sum((a.dtype.descr for a in arrays), []))
        if len(np.unique(newdtype.names)) != len(newdtype.names):
            raise ValueError('arrays have duplicate fields.')
        newrecarray = np.empty(len(arrays[0]), dtype=newdtype)
        for a in arrays:
            for name in a.dtype.names:
                newrecarray[name] = a[name]
        return newrecarray

    def model_weights(self, target_points, parnames=None, subinds=None,
                      itype='dt', force_triangulation=False, **extras):
        """Given an array of target coordinates and optionally a list
        of parnames, return the indices and weights of the model grid
        points that will interpolate to the target points.

        :param target_points: ndarray, shape(ntarg, npar)
            The points to which you want to interpolate.  Can be a
            numpy structured array of length NTARG or a simple numpy
            array of dimension [NTARG x NPAR], if parnames is also
            passed.  If a numpy structured array, then the field names
            should correspond to parameters in the `par` structure.
            1d arrays (i.e. a single target point) will be
            automatically upgraded to 2d arrays

        :param parnames: string sequence
            A list of the parameter names corresponding to the NPAR
            dimension of target_points.  These names should correspond
            to parameters in the `par` structure.

        :param subinds:
            An optional array identifying the members of the grid to
            be used for interpolation.

        :param itype: string (default 'dt')
            A string giving the type of interpolation to perform.
            Available options are:

            * 'dt' - use Delaunay triangulation.  This is suitable for
              irregular grids, though care must be taken that target
              points are within the convex hull described by the grid
            * 'idw' - inverse distance weighting to the NPAR nearest
              neighbors.
            * 'nearest' - nearest neighbor, as determined by a kd-tree
            * 'ndlinear' - not yet implemented'

        :param force_triangulation: bool (default True)
            Normally the triangulation and kdtree formed from the
            model grid is stored and reused unless the
            `triangle_dirtiness' attribute is greater than zero.  If
            this keyword is true, force the triangulation and kdtree
            to be regenerated regardless of 'dirtiness'

        :returns inds: ndarray, shape(ntarg, nind)
            An array that indexes the model parameter grid to give the
            interpolating models for each target point.  It has shape
            (NTARG, NIND) where nind is NPAR+1 for Delaunay
            triangulation, NPAR for 'idw', and 1 for 'nearest'

        :returns weights: ndarray, shape(ntarg, nind)
            The interpolation weight for each model grid point
            specified by inds, and of identical shape. The weights
            summed along the NIND direction will always be 1. Thus for
            'nearest' weights is always 1.
        """
        # deal with recarray input
        if parnames is None:
            parnames = target_points.dtype.names
            targets = np.array(target_points.tolist())
        else:
            targets = target_points
        targets = np.atleast_2d(targets)

        # If necessary rebuild the model delauynay triangulation and
        # KDTree.  Pull the grid points out of the model record data
        # and make an (nmodel,ndim) array of model grid parameter
        # values.  Need to loop to make sure order is correct.
        # if subinds is not None:
        #     force_triangulation = True
        # print( subinds)
        self.triangle_dirtiness += force_triangulation
        if self.triangle_dirtiness > 0:
            self.refresh_graphs(parnames, subinds=subinds)

        # pass the result to weightsDT
        if itype.lower() == 'dt':
            inds, weights = self.weightsDT(targets)
        elif itype.lower() == 'idw':
            inds, weights = self.weights_kNN_inverse_dist(targets, k=targets.shape[1])
        elif itype.lower() == 'nearest':
            inds, weights = self.weights_kNN_inverse_dist(targets, k=1)
        return inds, weights

    def refresh_graphs(self, parnames, subinds=None):
        """Given parameter names and optionally specific model
        indices, build a a Delaunay triangulation and a KDTree for the
        model points.
        """
        model_points = [np.squeeze(self.pars[subinds][pname])
                        for pname in parnames]
        model_points = np.array(model_points).transpose()  # (nmod, ndim)
        self.graphed_parameters = parnames
        self.triangle_dirtiness = 0
        # Delaunay triangulate
        self.dtri = Delaunay(model_points)
        # kdtree
        try:
            self.kdt = sklearn.neighbors.KDTree(model_points)
        except(NameError):
            pass

    def weightsDT(self, target_points):
        """The interpolation weights are determined from barycenter
        coordinates of the vertices of the enclosing Delaunay
        triangulation simplex. This allows for the use of irregular Nd
        grids. See also weights_1DLinear and
        weights_kNN_inverse_distance.

        :param target_points: ndarray, shape(ntarg,npar)
            The coordinates to which you wish to interpolate.

        :returns inds: ndarray, shape(ntarg,npar+1)
             The model indices of the interpolates.

        :returns weights: narray, shape (ntarg,npar+1)
             The weights of each model given by ind in the
             interpolates.
        """
        # Find the encompassing (hyper)triangle(s) for the desired
        # points, given a delauynay triangulation
        ndim = target_points.shape[-1]
        # Output triangle_inds is an (ntarg) array of simplex indices
        triangle_inds = self.dtri.find_simplex(target_points)
        # And get model indices of (hyper)triangle vertices. inds has
        # shape (ntarg,ndim+1)
        inds = self.dtri.vertices[triangle_inds, :]
        # get the barycenter coordinates through matrix multiplication
        # and dimensional juggling
        tmp = self.dtri.transform[triangle_inds, ndim, :]
        bary = np.dot(self.dtri.transform[triangle_inds, :ndim, :ndim],
                      (target_points - tmp).reshape(-1, ndim, 1))
        oned = np.arange(triangle_inds.shape[0])
        # ok.  in np 1.7 can add an axis to squeeze
        bary = np.atleast_2d(np.squeeze(bary[oned, :, oned, :]))
        # the last bary coordinate is 1-sum of the other coordinates
        last = 1 - bary.sum(axis=-1)
        weights = np.hstack((bary, last[:, np.newaxis]))

        outside = (triangle_inds == -1)
        weights[outside, :] = 0

        return inds, weights

        # loop implementation of the above for clarity
        #   npts = triangle_inds.shape[0]
        #   bary = np.zeros([npts,ndim+1])
        #   for i in xrange(npts):
        #      bary[i,:-1]= np.dot( dtri.transform[triangle_inds[0],:ndim,:ndim],
        #                         (target_points-dtri.transform[triangle_inds[i],ndim,:])

    def weights_kNN_inverse_dist(self, target_points, k=1):
        """The interpolation weights are determined from the inverse
        distance to the k nearest neighbors.

        :param target_points: ndarray, shape(ntarg,npar)
            The coordinates to which you wish to interpolate.

        :param k:
            The number of nearest neighbors to use.

        :returns inds: ndarray, shape(ntarg,npar+1)
             The model indices of the interpolates.

        :returns weights: narray, shape (ntarg,npar+1)
             The weights of each model given by ind in the
             interpolates.
        """
        dists, inds = self.kdt.query(target_points, k=k,
                                     return_distance=True)
        if k == 1:
            return inds, np.ones(inds.shape)
        weights = 1 / dists
        weights[np.isinf(weights)] = large_number
        weights = weights/weights.sum(axis=-1)
        return inds, weights

    def weights_1DLinear(self, model_points, target_points,
                         extrapolate=False, left=0.0, right=0.0):
        """The interpolation weights are determined from 1D linear
        interpolation.

        :param model_points: ndarray, shape(nmod)
            The parameter coordinate of the available models

        :param target_points: ndarray, shape(ntarg)
            The coordinate to which you wish to interpolate

        :returns inds: ndarray, shape(ntarg,2)
             The model indices of the interpolates

        :returns weights: narray, shape (ntarg,2)
             The weights of each model given by ind in the interpolates.
        """
        order = model_points.argsort()
        mod_sorted = model_points[order]

        x_new_indices = np.searchsorted(mod_sorted, target_points)
        x_new_indices = x_new_indices.clip(1, len(mod_sorted) - 1).astype(int)
        lo = x_new_indices - 1
        hi = x_new_indices
        x_lo = mod_sorted[lo]
        x_hi = mod_sorted[hi]
        width = x_hi - x_lo
        w_lo = (x_hi - target_points) / width
        w_hi = (target_points - x_lo) / width

        if extrapolate is False:
            # find places where target is above or below the model
            # range.
            above_scale = w_lo < 0
            below_scale = w_hi < 0
            # set the indices to be indentical in these cases
            lo[above_scale] = hi[above_scale]
            hi[below_scale] = lo[below_scale]
            # make the combined weights sum to ``left`` or ``right``
            w_lo[above_scale] = 0
            w_hi[above_scale] = left
            w_hi[below_scale] = 0
            w_lo[below_scale] = right

        inds = np.vstack([lo, hi]).T
        weights = np.vstack([w_lo, w_hi]).T

        return inds, weights

    def nearest_index(self, array, value):
        return (np.abs(array-value)).argmin(axis=-1)


class SpecLibrary(ModelLibrary):
    """Class to operate on spectral libraries. Methods are provided to
    interpolate the available model spectra (stored as a structured
    parameter array and a spectral array) to a certain set of
    parameters.  Subclasses are used to return the actual model
    spectrum given a set of model parameters.  Primary attributes are
    pars, spectra, wavelength.  Spectra should be of shape (NOBJ x
    NWAVE)
    """
    flux_unit = 'erg/s/cm^2/AA of 1solar mass at 10pc'

    def __init__(self):
        pass

    def spectra_from_pars(self):
        """This should take a numpy structured array of parameters of
        length nobj and return an (nobj, nwave) array of spectra
        """
        raise(NotImplementedError)

    def generateSEDs(self, pars, filterlist, wave_min=90, wave_max=1e7,
                     keepspec=False, intspec=False, attenuator=None, **extras):
        """
        :returns sed:
            ndarray of shape (nobj,nfilter)

        :returns lbol:
            ndarray of shape (nobj)

        :returns outspectra:
            ndarray of shape (nobj,nwave)
        """
        # don't use too much memory at once
        maxmod = 1e7/self.wavelength.shape[0]
        ngrid = pars.shape[0]
        sed = np.zeros([int(ngrid), len(filterlist)])
        lbol = np.zeros(int(ngrid))
        outspectra = None
        if keepspec:
            outspectra = np.zeros([ngrid, self.wavelength.shape[0]])
        elif intspec:
            outspectra = np.zeros(self.wavelength.shape[0])

        # split big model grids to avoid memory constraints
        i = 0
        while (i*maxmod <= ngrid):
            s1, s2 = int((i)*maxmod), int(np.min([(i+1)*maxmod-1, ngrid]))
            spec = self.spectra_from_pars(pars[int(s1):int(s2)], **extras)
            if attenuator is not None:
                spec = attenuator.attenuate_spectrum(self.wavelength, spec,
                                                     pars[s1:s2], **extras)
            sed[s1:s2, :] = observate.getSED(self.wavelength, spec, filterlist)
            lbol[s1:s2] = observate.Lbol(self.wavelength, spec, wave_min=wave_min,
                                         wave_max=wave_max)
            i += 1
            if keepspec:
                outspectra[s1:s2, :] = spec
            elif intspec:
                outspectra += spec.sum(axis=0)

        return sed, lbol, outspectra

    def interpolate_to_pars(self, target_points, parnames=None,
                            subinds=None, itype='dt', **extras):
        """Method to obtain the model spectrum for a given set of
        parameter values via interpolation of the model grid. The
        interpolation weights are determined from barycenters of a
        Delaunay triangulation or nLinear interpolation or
        k-nearest-neighbor inverse distance.

        :param target_points: ndarray, shape (ntarg,ndim)
            Desired model parameters. Can also be a structured array
            with fields named for the model parameters.  See
            ModelLibrary.model_weights()

        :param subinds: ndarray
            Indices of the model pars structure to use in
            interpolation. Allows for only portions of the model
            library to be used.

        :param parnames: string sequence, len (ndim)
            The names of the model library parameters

        :returns spectra: ndarray, shape (ntarg, nwave)
            The interpolated spectra.
        """
        inds, weights = self.model_weights(target_points, parnames=parnames,
                                           itype=itype, subinds=subinds, **extras)
        if subinds is not None:
            inds = subinds[inds]
        return self.combine_weighted_spectra(inds, weights)

    def combine_weighted_spectra(self, inds, weights):
        """Weight self.spectra using broadcasting, then sum the
        weighted spectra.  Should add a switch to sum the weights
        first so that the output is a single spectrum.

        :param inds: shape (ntarg,nint)
            Indices of the models to sum

        :param weights: shape (ntarg, nint)
            Weights of the models corresponding to inds.

        :returns spec: shape (nwave,ntarg).
            The weighted sum.
        """
        tmp = weights * (self.spectra[inds].transpose(2, 0, 1))
        return (tmp.sum(axis=2)).T

    def read_model_from_fitsbinary(self, filename, parnames,
                                   wavename='WAVE', fluxname='F_LAMBDA'):
        """Read spectra from a fits binary table in a certain format.
        """
        # if os.ispath(filename) is False:
        #      raise IOError('read_model_from_fitsbinary:
        #                    ',filename,' does not exist')
        fits = pyfits.open(filename)
        # Parse the FITS recarray and assign ModelGrid parameter,
        # spectra, and wavelength attributes
        wavelength = fits[1].data[0][wavename]
        spectra = fits[1].data[fluxname]  # (nmod,nwave)
        # Need to loop over pars and restruct to get the desired
        # order.  Really?  Not really.  order is unimportant
        pars, partype = [], []
        for pname in parnames:
            pars.append(np.squeeze(fits[1].data[pname]))
            partype.append(fits[1].data[pname].dtype)  # fix
        # list ->(nmod, npar) -> Structured array
        pars = self.structure_array(np.array(pars).transpose(),
                                    parnames, types=partype)
        fits.close()
        return wavelength, spectra, pars
