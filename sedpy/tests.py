import numpy as np

from . import observate


def _get_all_filters():
    fnames = np.array(observate.list_available_filters())
    allfilters = observate.load_filters(fnames)
    return allfilters


def test_filter_load():
    allfilters = _get_all_filters()
    w = [f.wave_effective for f in allfilters]
    
    
def test_gridded_filters():
    allfilters = np.array(_get_all_filters())
    w = np.array([f.wave_effective for f in allfilters])
    fnames = np.array([f.name for f in allfilters])
    
    good = w < 1e5
    obs = {}
    obs['filters'] = allfilters[good][0:40]

    spec = np.random.uniform(0, 1.0, 5996)
    wave = np.exp(np.linspace(np.log(90), np.log(1e6), len(spec)))
    #%timeit m = observate.getSED(wave, spec, obs['filters'])
    #3.22 ms
    m_default = observate.getSED(wave, spec, obs['filters'])

    wlo, whi, dlo = [], [], []
    for f in obs['filters']:
        dlnlam = np.gradient(f.wavelength)/f.wavelength
        wlo.append(f.wavelength.min())
        dlo.append(dlnlam.min())
        whi.append(f.wavelength.max())
    wmin = np.min(wlo)
    wmax = np.max(whi)
    dlnlam = np.min(dlo)
    obs['filters'] = observate.load_filters(fnames[good][0:40],
                                            dlnlam=dlnlam, wmin=wmin)

    lnlam = np.exp(np.arange(np.log(wmin), np.log(wmax), dlnlam))
    lnspec = np.interp(lnlam, wave, spec)

    m_grid = m = observate.getSED(lnlam, lnspec, obs['filters'],
                                  gridded=True)
    # make sure good to 5% (most far better, hipparcos_B is a problem)
    assert np.allclose(m_grid, m_default, atol=5e-2)


def test_filter_properties():
    """Compare to the values obtained from the K-correct code
    (which uses a slightly different Vega spectrum)
    """
    filternames = ['galex_FUV',
                   'sdss_u0','sdss_g0','sdss_r0','sdss_i0',
                   'spitzer_irac_ch2']
    weff_kcorr = [1528.0,
                  3546.0, 4669.6, 6156.25, 7471.57,
                  44826.]
    msun_kcorr = [18.8462,
                  6.38989, 5.12388, 4.64505, 4.53257,
                  6.56205]
    ab2vega_kcorr = [2.3457,
                     0.932765, -0.0857, 0.155485, 0.369598,
                     3.2687]

    filterlist = observate.load_filters(filternames)
    for i, f in enumerate(filterlist):
        #print(filterlist[i].wave_effective, filterlist[i].solar_ab_mag,
        #      filterlist[i].ab_to_vega)
        assert (abs(f.wave_effective - weff_kcorr[i]) <
                (weff_kcorr[i] * 0.01))
        assert abs(f.solar_ab_mag - msun_kcorr[i]) < 0.05
        # the below fails because of the vega spectrum used by k_correct
        # assert abs(f.ab_to_vega+ab2vega_kcorr[i]) < 0.05
