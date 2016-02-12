Frequently Asked Questions
============

How do I add filter transmission curves?
--------
Each transmission curve should be in its own file named ``filter_name.par``,
where ``filter_name`` is a unique identifier that will be used to identify the filter in the code
(e.g., ``galex_FUV.par``).

There are two options for the file format and location, one simple and one preferred and more stable.

*Simple:*

Make ``filter_name.par`` a two-column space delimited ascii file where
the first column is wavelength in angstroms and
the second column is detector signal per photon (arbitrarily normalized).
Comments must be preceded by a ``#``.
Then you can access it in python easily with e.g.

.. code-block:: python

		from sedpy import observate
		filt = observate.Filter("filter_name", directory="/path/to/the/file/")
		filt_list = observate.load_filters(["filter_name"], directory="/path/to/the/file/")
		
If using ``load_filters`` then all transmission files must be in the same directory.
Optionally you can place the files you create  in the ``<download_dir>/sedpy/sedpy/data/filters/`` directory
and reinstall with ``cd <download_dir>/sedpy; python setup.py install``.
Then you don't need to specify the ``directory`` keyword.

*Preferred:*

The preferred format of the filter transmission file was inherited from the yanny files used
for `Kcorrect <http://howdy.physics.nyu.edu/index.php/Kcorrect>`_ and is:

.. code-block:: C
		
		# comments
		# more comments

		typedef struct {
		  double lambda;
		  double pass;
		} KFILTER;

		KFILTER    1000.0   0.0
		KFILTER    1001.0   1.1540
		.
		.
		.
		KFILTER    3000.0   0.0
		
where the second column is the wavelength in Angstroms and the third column is detector signal per photon.
   
Once you have a correctly formatted transmission file,
simply place this file in the ``<download_dir>/sedpy/sedpy/data/filters/`` directory and reinstall with ``cd <download_dir>/sedpy; python setup.py install``

The new filter will now be accesible using ``filter_name``

.. code-block:: python

		from sedpy import observate
		filt = observate.Filter("filter_name")
