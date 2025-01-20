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

For some filters, multiple versions of the transmission (e.g. corresponding to
different chips) are stored in the same file in different columns.  While the
default column is always the first one, other columns can be accessed by name as follows

```python
from sedpy import observate
# filter corresponding to an average over SCAs
filt_mean =  observate.Filter("roman_wfi_f146")
# filter corresponding to SCA1
filt_sca01 = observate.Filter("roman_wfi_f146", trans_colname="sca01")

print(filt_mean.name, filt_mean.wave_effective)
print(filt_sca01.name, filt_sca01.wave_effective)

```


Spectral smoothing
------------------
