sedpy
=====

Modules for storing and operating on astronomical source spectral energy distributions.

Installation & setup:
---------------------

.. code-block:: shell

		git clone https://github.com/bd-j/sedpy
		cd sedpy
		pip install .

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
------------

* ``observate`` has methods for generating synthetic photometry through any filters,
  and classes for dealing with filters generally. There is some functionality for spectra
  (vaccum to air conversions).
  With a huge debt to Mike Blanton's `kcorrect <https://github.com/blanton144/kcorrect>`_ code .

* ``attenuation`` contains simple dust attenuation methods.

* ``extinction`` (Deprecated) classes for a detailed modeling of extinction curves,
  following the Fitzpatrick & Massa parameterizations.
  See `dust_extinction <https://dust-extinction.readthedocs.io/en/stable/>`_ instead.

* ``photometer`` (Deprecated) has some basic aperture photometry algorithms.
  See `photutils <https://photutils.readthedocs.io/en/stable/>`_ instead.

* ``ds9region`` (Deprecated) has some simple ds9 region classes.

* ``modelgrid`` (Deprecated) is a module with classes for the storage and interpolation of
  model SEDs. Largely superceded by ``scipy.interpolate`` algorithms.
