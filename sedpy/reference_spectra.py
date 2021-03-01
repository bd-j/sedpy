# -*- coding: utf-8 -*-

"""reference_spectra.py  - load and store Vega and Solar reference spectra.

Units of the reference spectra are AA and erg/s/cm^2/AA
"""

import numpy as np
import os
try:
    from pkg_resources import resource_filename
except ImportError:
    pass
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits

__all__ = ["vega", "solar", "sedpydir"]

# ----------
# Load useful reference spectra
# ----------

sedpydir, f = os.path.split(__file__)
sedpydir = sedpydir
try:
    vega_file = resource_filename('sedpy', 'data/alpha_lyr_stis_005.fits')
except:
    vega_file = os.path.join(sedpydir, 'data', 'alpha_lyr_stis_005.fits')

# This file should be in AA and erg/s/cm^2/AA
if os.path.isfile(vega_file):
    fits = pyfits.open(vega_file)
    vega = np.column_stack((fits[1].data.field('WAVELENGTH'),
                            fits[1].data.field('FLUX')))
    fits.close()
else:
    raise ValueError('Could not find Vega '
                     'spectrum at {0}'.format(vega_file))

try:
    solar_file = resource_filename('sedpy', 'data/sun_kurucz93.fits')
except:
    solar_file = os.path.join(sedpydir, 'data', 'sun_kurucz93.fits')

# conversion to d=10 pc from 1 AU
rat = (1.0 / (3600 * 180 / np.pi * 10))**2.0

# This file should be in AA and erg/s/cm^2/AA at 1AU
if os.path.isfile(solar_file):
    fits = pyfits.open(solar_file)
    solar = np.column_stack((fits[1].data.field('WAVELENGTH'),
                             fits[1].data.field('FLUX') * rat))
    fits.close()
else:
    raise ValueError('Could not find Solar '
                     'spectrum at {0}'.format(solar_file))
