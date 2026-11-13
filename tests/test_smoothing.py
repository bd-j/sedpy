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
