sedpy
======

Modules for storing and operating on astronomical source spectral energy distributions.

Installation & setup:
--------------

.. code-block:: shell

		git clone https://github.com/bd-j/sedpy
		cp /path/to/your/favorite/filters/*par sedpy/sedpy/data/filters/
		cd sedpy
		python setup.py install

Then in python, e.g.,

.. code-block:: python

		from sedpy import observate
		# get magnitude from a spectrum:
		filt = observate.Filter("sdss_r0")
		mag = filt.ab_mag(angstroms, f_lambda_cgs)
		# or get several magnitudes at once
		filterlist = observate.load_filters(["galex_NUV", "sdss_r0"])
		mags = observate.getSED(angstroms, f_lambda_cgs, filterlist=filters)

see the `FAQ`_

.. _FAQ: docs/faq.rst


Description:
---------------

* ``observate`` has methods for generating synthetic photometry through any filters,
  and classes for dealing with filters generally. There is some functionality for spectra (vaccum to air conversions).
  With a huge debt to Mike Blanton's `kcorrect <https://github.com/blanton144/kcorrect>`_ code .

* ``attenuation`` contains simple dust attenuation methods.

* ``extinction`` contains classes for a detailed modeling of extinction curves,
  following the Fitzpatrick and Massa parameterizations.
  (Deprecated, see `dust_extinction <https://dust-extinction.readthedocs.io/en/stable/>`_)

* ``photometer`` has some basic aperture photometry algorithms.
   (Deprecated, see `photutils <https://photutils.readthedocs.io/en/stable/>`_)

* ``ds9region`` has some simple ds9 region classes. (Deprecated)

* ``modelgrid`` is a module with classes for the storage and interpolation of
  model SEDs.  (Deprecated, largely superceded by ``scipy.interpolate`` algorithms)

* ``yanny`` (from Erin Sheldon) is used internally for reading filter curves,
  included here for convenience.
