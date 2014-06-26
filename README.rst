sedpy
======

Modules for storing and operating on astronomical source spectral energy distributions.

Installation & setup:

1. ``git clone http://blah``
2. ``cd sedpy; python setup.py install``
3. cp any desired filters from the k_correct filter directory to
   sedpy/data/filters/
4. then in python, e.g., ``from sedpy import observate``
    
Description:

* ``observate`` has methods for generating synthetic photometry
  through any filters, and classes for dealing with filters
  generally. Functionality for spectra is being added slowly. With a
  huge debt to Mike Blanton's k_correct code .
  
* ``attenuation`` contains dust attenuation methods (very
  preliminary).

* ``extinction`` contains classes for a detailed modeling of
  extinction curves, following the Fitzpatrick and Massa
  parameterizations.

* ``modelgrid`` is a module with classes for the storage and
  interpolation of model SEDs (using numpy structured arrays, usually,
  and Delaunay triangulation or inverse distance weighting of the
  k-nearest neighbors and KD-trees)

* ``photometer`` has some basic aperture photometry algorithms.

* ``ds9region`` has some simple ds9 region classes.

* ``yanny`` (from Erin Sheldon) is for reading filter curves, included here for convenience
