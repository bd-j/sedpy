sedpy
======

Modules for storing and operating on astronomical source spectral energy distributions.

Installation & setup:

1. ``git clone https://github.com/bd-j/sedpy``
   
2. copy any desired filters from the ``kcorrect`` filter directory (or your
   own filter transmission curves) to ``sedpy/sedpy/data/filters/``
   
3. ``cd sedpy; python setup.py install``
   
4. then in python, e.g., ``from sedpy import observate``

see the [FAQ](docs/faq.rst)
   
Description:

* ``observate`` has methods for generating synthetic photometry through any filters,
  and classes for dealing with filters generally.
 There is some functionality for spectra (vaccum to air conversions).
  With a huge debt to Mike Blanton's `kcorrect <https://github.com/blanton144/kcorrect>`_ code .
  
* ``attenuation`` contains simple dust attenuation methods.

* ``extinction`` contains classes for a detailed modeling of extinction curves,
  following the Fitzpatrick and Massa parameterizations.

* ``modelgrid`` is a module with classes for the storage and interpolation of
  model SEDs (using numpy structured arrays, usually, and Delaunay triangulation or
  inverse distance weighting of the k-nearest neighbors and KD-trees)

* ``photometer`` has some basic aperture photometry algorithms.

* ``ds9region`` has some simple ds9 region classes.

* ``yanny`` (from Erin Sheldon) is for reading filter curves, included here for convenience

* ``smoothing.py`` contains obsolete methods for spectral smoothing.  See
  `prospect.utils.smoothing <https://github.com/bd-j/prospector/blob/master/prospect/utils/smoothing.py>`_ instead.
