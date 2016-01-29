Frequently Asked Questions
============

How do I add filter transmission curves?
--------
First you must put your transmission curves in the proper file format.
Each transmission curve should be in its own file named ``filter_name.par``,
where ``filter_name`` is a unique identifier that will be used to identify the filter in the code
(e.g., ``galex_FUV.par``).
It should not have any spaces.

The format of the filter transmission file was inherited from the yanny files used
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
