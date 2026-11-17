"""Failures in the band-analysis numerics must be visible, not silently masked."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("matplotlib")

from rsutil.bands.analysis.magnetic.effective_mass import EffectiveMassCalculator
from rsutil.bands.analysis.topological.berry_curvature import BerryCurvatureCalculator
from rsutil.bands.dos.mesh import FullBZCalculator
from rsutil.bands.readers import detect_hamiltonian_format, detect_lattice_format


def _calculator() -> EffectiveMassCalculator:
    # Bypass __init__ (extremum search): only the fitting helper is exercised.
    return EffectiveMassCalculator.__new__(EffectiveMassCalculator)


def test_quadratic_fit_success_has_no_warning(recwarn):
    k = np.linspace(-0.2, 0.2, 9)
    e0, curv, r2 = _calculator()._quadratic_fit(k, 1.5 + 3.0 * k**2)
    assert e0 == pytest.approx(1.5) and curv == pytest.approx(3.0) and r2 == pytest.approx(1.0)
    assert not [w for w in recwarn if issubclass(w.category, RuntimeWarning)]


def test_quadratic_fit_failure_warns_and_returns_nan():
    # NaN input makes scipy raise ValueError: previously swallowed silently.
    k = np.linspace(-0.2, 0.2, 9)
    with pytest.warns(RuntimeWarning, match="Quadratic band fit failed"):
        e0, curv, r2 = _calculator()._quadratic_fit(k, np.full_like(k, np.nan))
    assert np.isnan(e0) and np.isnan(curv) and r2 == 0.0


def test_unexpected_exceptions_are_not_swallowed():
    # Programming errors (and Ctrl-C) must propagate, which a bare except hid.
    calc = _calculator()

    def boom(*_args, **_kwargs):
        raise KeyboardInterrupt

    calc.calculate_effective_mass_at_extremum = boom
    calc.extrema = [1]
    with pytest.raises(KeyboardInterrupt):
        calc.calculate_all_effective_masses()


def test_all_effective_masses_skips_failed_extremum_with_warning():
    calc = _calculator()
    calc.extrema = [1, 2]

    def fake(idx, _k_range):
        if idx == 1:
            raise ValueError("degenerate")
        return {"goodness_of_fit": 0.99, "effective_mass": 0.5}

    calc.calculate_effective_mass_at_extremum = fake
    with pytest.warns(RuntimeWarning, match="Skipping extremum 1"):
        results = calc.calculate_all_effective_masses()
    assert len(results) == 1 and results[0]["effective_mass"] == 0.5


def test_chern_histogram_marks_failed_band_as_missing_not_zero(tmp_path):
    calc = BerryCurvatureCalculator.__new__(BerryCurvatureCalculator)
    calc.nbands = 3

    def chern(band):
        if band == 1:
            raise ValueError("singular")
        return 1.0 if band == 0 else 0.0

    calc.calculate_chern_number = chern
    out = tmp_path / "chern.png"
    with pytest.warns(RuntimeWarning, match="band 1 could not be computed"):
        calc.plot_chern_number_histogram(str(out))
    assert out.is_file()


def test_tetrahedra_fallback_is_announced(monkeypatch):
    import scipy.spatial

    def broken(*_a, **_k):
        raise ValueError("qhull failed")

    monkeypatch.setattr(scipy.spatial, "Delaunay", broken)
    mesh = FullBZCalculator.__new__(FullBZCalculator)
    with pytest.warns(RuntimeWarning, match="simplified tetrahedra"):
        tetra = mesh._generate_tetrahedra(np.random.default_rng(0).random((8, 3)), (2, 2, 2))
    assert len(tetra) > 0


def test_format_detection_survives_unreadable_files(tmp_path):
    binary = tmp_path / "blob.dat"
    binary.write_bytes(b"\xff\xfe\x00\x80" * 64)
    assert detect_lattice_format(str(binary)) == "rs-lmto"
    assert detect_hamiltonian_format(str(binary)) == "rs-lmto"
    assert detect_lattice_format(str(tmp_path / "missing.dat")) == "rs-lmto"
