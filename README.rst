sedpy
=====


Modules for storing and operating on astronomical source spectral energy distributions.

.. image:: https://github.com/bd-j/sedpy/workflows/Tests/badge.svg
  :target: https://github.com/bd-j/sedpy/actions?query=workflow%3ATests

.. image:: https://readthedocs.org/projects/sedpy/badge/?version=latest
    :target: https://sedpy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Installation & setup:
---------------------
``sedpy`` is pip installable:

.. code-block:: shell

		python -m pip install astro-sedpy

Or you can install the latest version from github:

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

For the filters available by default see the `filter_list`_.
For adding transmission curves, see these `docs`_.

.. _filter_list: sedpy/data/filters/README.md
.. _docs: docs/transmissions.rst

This code can be referenced as:

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4582723.svg
   :target: https://doi.org/10.5281/zenodo.4582723

Description:
------------

* ``observate`` has methods for generating synthetic photometry through any filters,
  and classes for dealing with filters generally. There is some functionality for spectra
  (vaccum to air conversions).
  With a huge debt to Mike Blanton's `kcorrect <https://github.com/blanton144/kcorrect>`_ code .

* ``attenuation`` contains simple dust attenuation methods.

* ``smoothing`` methods for smoothing well-sampled spectra.

* ``extinction`` (Deprecated) classes for a detailed modeling of extinction curves,
  following the Fitzpatrick & Massa parameterizations.
  See `dust_extinction <https://dust-extinction.readthedocs.io/en/stable/>`_ instead.
