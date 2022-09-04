Quickstart
==========

Filter projections
------------------

For a given spectrum with matched arrays of wavelength (given in Angstrom) and
flux (given in erg/s/cm^2/s) you can do as follows for a given filter:

```python
from sedpy import observate
# get AB magnitude from a spectrum:
filt = observate.Filter("sdss_r0")
mag = filt.ab_mag(angstroms, f_lambda_cgs)
# or get several magnitudes at once
filterlist = observate.load_filters(["galex_NUV", "sdss_r0"])
mags = observate.getSED(angstroms, f_lambda_cgs, filterlist=filters)
```

Spectral smoothing
------------------

