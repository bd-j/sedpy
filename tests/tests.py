# -*- coding: utf-8 -*-

import numpy as np
from numpy.testing import assert_allclose

from sedpy import observate


def _get_all_filters():
    fnames = observate.list_available_filters()
    allfilters = observate.load_filters(fnames)
    return allfilters


def test_filter_load():
    fnames = observate.list_available_filters()
    for fn in fnames:
        try:
            f = observate.Filter(fn)
        except(IOError, ValueError, IndexError, AssertionError) as e:
            print(f"Could not build `{fn}`")
            raise(e)
        w = f.wave_effective
        assert np.isfinite(w), f"effective wavelength undefined for {fn}"


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

    m_grid = observate.getSED(lnlam, lnspec, obs['filters'],
                              gridded=True)
    # make sure good to 5% (most far better, hipparcos_B is a problem)
    assert np.allclose(m_grid, m_default, atol=5e-2)


def test_filterset():
    from sedpy.observate import FilterSet
    flist = [f for f in observate.list_available_filters()
             if ("jw" in f) or ("wise" in f) or ("galex" in f)]
    filterset = FilterSet(flist, dlnlam=1e-4)
    source = observate.vega.T

    mags = observate.getSED(source[0], source[1], filterset)
    omags = observate.getSED(source[0], source[1], observate.load_filters(flist))

    good = np.isclose(mags, omags, atol=1.5e-3)
    assert np.all(good), f"filters {np.array(flist)[~good]} failed with {(mags-omags)[~good]}"


def test_shapes():

    from sedpy.observate import FilterSet
    flist = [f for f in observate.list_available_filters()
             if ("jw" in f) or ("wise" in f) or ("galex" in f)]
    sw, sf0 = observate.vega.T
    assert sf0.ndim == 1

    for fl in [flist[:1], flist]:
        Nf = len(fl)
        filterset = FilterSet(fl)

        # test input flux of shape (nwave,)
        mags = observate.getSED(sw, sf0, filterset)
        assert mags.shape == (Nf,), "getSED with FilterSet gives {}, expected {}".format(mags.shape, (Nf,))
        omags = observate.getSED(sw, sf0, observate.load_filters(fl))
        assert omags.shape == (Nf,), "getSED with Filter list gives {}, expected {}".format(omags.shape, (Nf,))

        # test vector inputs, though inputs of shape (1, nwave) get turned into (nfilt,)
        for N in [1, 10]:
            if N == 1:
                truth = (Nf,)
            else:
                truth = (N, Nf)

            sf = np.tile(sf0, (N, 1))
            mags = observate.getSED(sw, sf, filterset)
            assert mags.shape == truth, "getSED with FilterSet gives {}, expected {}".format(mags.shape, truth)
            omags = observate.getSED(sw, sf, observate.load_filters(fl))
            assert omags.shape == truth, "getSED with Filter list gives {}, expected {}".format(omags.shape, truth)

            assert np.allclose(mags, omags, 5e-3), mags - omags


def test_gridded_shapes():

    from sedpy.observate import FilterSet

    # array of object spectra
    Nobj, Nwave = 100, 3000
    spec = np.ones([Nobj, Nwave], dtype=float)
    wave = np.exp(np.linspace(np.log(0.5e4), np.log(5.1e4), Nwave))

    fnames = [f"jwst_f{b}" for b in ["070w", "090w", "115w", "150w", "200w", "335m"]]
    fnames.sort()

    wmin = wave.min()
    dlnlam = np.gradient(np.log(wave))
    dlnlam_filters = dlnlam.min()
    assert np.allclose(dlnlam, dlnlam[0])

    filterset = FilterSet(fnames, dlnlam=dlnlam_filters, wmin=wmin)
    inds = np.log(wave) <= (np.log(filterset.lam.max()) + dlnlam_filters*0.1)

    maggies = filterset.get_sed_maggies(spec[:, inds])
    assert maggies.shape == (100, len(fnames))


def test_filter_cols():
    from sedpy import observate
    # filter corresponding to an average over SCAs
    filt_mean =  observate.Filter("roman_wfi_f146")
    # filter corresponding to SCA1
    filt_sca01 = observate.Filter("roman_wfi_f146", trans_colname="sca01")
    k = filt_mean.name, filt_mean.wave_effective
    j = filt_sca01.name, filt_sca01.wave_effective
    assert (filt_mean.name != filt_sca01.name)


def test_filter_properties():
    """Compare to the values obtained from the K-correct code
    (which uses a slightly different Vega spectrum)
    """
    filternames = ['galex_FUV',
                   'sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0',
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
