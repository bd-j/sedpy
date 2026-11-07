# -*- coding: utf-8 -*-

import numpy as np
from numpy.testing import assert_allclose

from sedpy.smoothing import smooth_fft, smoothspec

CKMS = 2.998e5


def _log_grid(lo=4700.0, hi=5300.0, dv=10.0):
    """A grid uniform in log wavelength with steps of ``dv`` km/s."""
    return np.exp(np.arange(np.log(lo), np.log(hi), dv / CKMS))


def _absorption_spectrum(wave, width=12.7, nline=60, seed=0):
    """Unit continuum with Gaussian absorption lines of ``width`` km/s, barely sampled."""
    rng = np.random.default_rng(seed)
    v = CKMS * np.log(wave[:, None] / rng.uniform(4800.0, 5200.0, nline)[None, :])
    return 1.0 - np.sum(rng.uniform(0.05, 0.6, nline) * np.exp(-0.5 * (v / width) ** 2), axis=1)


def test_smooth_fft_keeps_length():
    """Odd and even inputs come back with their own length."""
    for n in (1024, 1023):
        spec = np.random.default_rng(1).normal(size=n)
        assert smooth_fft(1.0, spec, 3.0).size == n


def test_vel_fft_matches_direct_on_log_uniform_grid():
    """On a log-uniform grid the FFT path matches the direct convolution of smooth_vel."""
    wave = _log_grid()
    spec = _absorption_spectrum(wave)
    core = (wave > 4850.0) & (wave < 5150.0)
    for sigma in (30.0, 150.0):
        direct = smoothspec(wave, spec, sigma, outwave=wave, smoothtype="vel", fftsmooth=False)
        fft = smoothspec(wave, spec, sigma, outwave=wave, smoothtype="vel", fftsmooth=True)
        assert_allclose(fft[core], direct[core], rtol=0, atol=1e-6)


def test_vel_fft_odd_length_log_uniform_grid():
    """An odd number of log-uniform samples is smoothed and returned at full length."""
    wave = _log_grid()
    wave = wave[: wave.size - (1 - wave.size % 2)]
    assert wave.size % 2 == 1
    out = smoothspec(wave, _absorption_spectrum(wave), 60.0, outwave=wave, smoothtype="vel")
    assert out.size == wave.size and np.all(np.isfinite(out))


def test_vel_fft_non_uniform_grid_is_resampled():
    """A grid that is not log-uniform still goes through the resampled FFT path."""
    rng = np.random.default_rng(2)
    wave = np.sort(rng.uniform(4700.0, 5300.0, 4000))
    spec = 1.0 + 0.1 * np.sin(wave / 3.0)
    out = smoothspec(wave, spec, 60.0, outwave=wave, smoothtype="vel", fftsmooth=True)
    direct = smoothspec(wave, spec, 60.0, outwave=wave, smoothtype="vel", fftsmooth=False)
    core = (wave > 4850.0) & (wave < 5150.0)
    assert out.size == wave.size
    # the resampled path interpolates first, so it agrees only to a few 1e-3
    assert_allclose(out[core], direct[core], rtol=0, atol=5e-3)


def test_smooth_fft_upsample_keeps_samples():
    """upsample=1 equals smooth_fft; with upsample=4 every 4th point is a smoothed input sample."""
    from sedpy.smoothing import smooth_fft_upsample

    for n in (1024, 1023):
        spec = np.random.default_rng(3).normal(size=n)
        base = smooth_fft(10.0, spec, 20.0)
        assert_allclose(smooth_fft_upsample(10.0, spec, 20.0, 1), base, rtol=0, atol=1e-14)
        fine = smooth_fft_upsample(10.0, spec, 20.0, 4)
        assert fine.size == 4 * (n - 1) + 1
        assert_allclose(fine[::4], base, rtol=0, atol=1e-12)


def _analytic(wave, sigma, width=12.7, nline=60, seed=0):
    """The spectrum of _absorption_spectrum smoothed analytically by ``sigma`` km/s."""
    rng = np.random.default_rng(seed)
    total = np.hypot(width, sigma)
    v = CKMS * np.log(wave[:, None] / rng.uniform(4800.0, 5200.0, nline)[None, :])
    depth = rng.uniform(0.05, 0.6, nline) * width / total
    return 1.0 - np.sum(depth * np.exp(-0.5 * (v / total) ** 2), axis=1)


def test_vel_fft_fine_outwave_is_upsampled():
    """On a finer outwave, Fourier upsampling keeps narrow lines; linear interpolation does not.

    Lines of 12.7 km/s on 10 km/s pixels, smoothed by 13-25 km/s (a stellar LOSVD of 18-28 km/s
    with the 12.7 km/s library resolution removed), evaluated on a 0.05 A grid.
    """
    wave = _log_grid(4500.0, 5500.0)
    spec = _absorption_spectrum(wave)
    fine = np.arange(4850.0, 5150.0, 0.05)
    for losvd in (18.0, 23.0, 28.0):
        sigma = np.sqrt(losvd**2 - 12.7**2)
        truth = _analytic(fine, sigma)
        out = smoothspec(wave, spec, sigma, outwave=fine, smoothtype="vel")
        native = smoothspec(wave, spec, sigma, outwave=wave, smoothtype="vel")
        linear = np.interp(fine, wave, native)
        assert np.max(np.abs(out - truth)) < 2e-3
        assert np.max(np.abs(linear - truth)) > 5e-3
