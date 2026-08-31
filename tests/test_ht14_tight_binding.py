"""HT-14 versioned TB parsing, reconstruction and deterministic ML export."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from rsutil.bands.readers import HamiltonianReader, LatticeReader, load_hamiltonian
from rsutil.results import (
    TightBindingModel,
    export_ml_dataset,
    read_fable_v4_tb,
    tight_binding_from_legacy_reader,
    tight_binding_from_ml_dataset,
)
from rsutil.results.tight_binding import _match_legacy_translation
from rsutil.core import RSStructure, StructureSite


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "tb" / "fable_v4" / "rs_tb"


def test_fable_v4_versioned_pair_preserves_identity_and_reconstructs_hk() -> None:
    model = read_fable_v4_tb(FIXTURE.with_suffix(".meta"))

    assert isinstance(model, TightBindingModel)
    assert model.energy_unit.value == "eV"
    assert model.spin_representation == "spinor"
    reverse = next(term for term in model.hoppings if term.R == (-1, 0, 0))
    assert (reverse.site_i, reverse.site_j, reverse.value) == (1, 0, 0.2 - 0.1j)

    H = model.hamiltonian((0.25, 0.0, 0.0))
    expected = np.zeros((4, 4), dtype=complex)
    expected[0, 0] = 0.1
    expected[1, 1] = 0.2
    expected[2, 2] = -0.1
    expected[3, 3] = 0.3
    expected[0, 2] = (0.2 + 0.1j) * np.exp(0.5j * np.pi)
    expected[2, 0] = (0.2 - 0.1j) * np.exp(-0.5j * np.pi)
    np.testing.assert_allclose(H, expected)
    np.testing.assert_allclose(model.eigenvalues((0.25, 0.0, 0.0)), np.linalg.eigvalsh(expected))


def test_ml_export_is_complex_lossless_and_byte_deterministic(tmp_path: Path) -> None:
    model = read_fable_v4_tb(FIXTURE.with_suffix(".meta"))
    first_npz, first_json = export_ml_dataset(model, tmp_path / "first")
    second_npz, second_json = export_ml_dataset(model, tmp_path / "second")

    assert hashlib.sha256(first_npz.read_bytes()).hexdigest() == hashlib.sha256(second_npz.read_bytes()).hexdigest()
    assert first_json.read_text(encoding="utf-8") == second_json.read_text(encoding="utf-8")
    restored = tight_binding_from_ml_dataset(first_npz)
    np.testing.assert_allclose(restored.hamiltonian((0.0, 0.0, 0.0)), model.hamiltonian((0.0, 0.0, 0.0)))
    assert restored.hoppings == model.hoppings


def test_legacy_reader_remains_available_and_conversion_requires_lattice_identity() -> None:
    lattice_path = ROOT / "examples" / "bands" / "FeCo" / "lattice.nml"
    ham_path = ROOT / "examples" / "bands" / "FeCo" / "rs2txt.dat"
    legacy = load_hamiltonian(str(ham_path))
    assert isinstance(legacy, HamiltonianReader)
    model = tight_binding_from_legacy_reader(legacy, LatticeReader(str(lattice_path)))
    assert model.energy_unit.value == "unknown"
    assert model.dimension == 36
    assert model.hoppings


def test_legacy_translation_solver_accepts_long_range_periodic_blocks() -> None:
    structure = RSStructure(
        lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        length_unit="angstrom",
        coordinate_convention="cartesian",
        sites=(StructureSite(0, (0, 0, 0), "X", None),),
    )

    assert _match_legacy_translation(structure, 0, (-6.0, 0.0, 0.0)) == (6, 0, 0, 0)


def test_control_nsp_records_soc_and_magnetic_representation() -> None:
    model_nsp1 = read_fable_v4_tb(
        FIXTURE.with_suffix(".meta"),
        input_path=ROOT / "examples" / "B2Fe" / "input.nml",
    )
    model_nsp2 = read_fable_v4_tb(
        FIXTURE.with_suffix(".meta"),
        input_path=ROOT / "examples" / "FePt" / "input.nml",
    )
    model_nsp4 = read_fable_v4_tb(
        FIXTURE.with_suffix(".meta"),
        input_path=ROOT / "examples" / "Mn3Sn" / "input.nml",
    )

    assert (model_nsp1.control_nsp, model_nsp1.soc_enabled, model_nsp1.magnetic_representation) == (
        1,
        False,
        "collinear_magnetic",
    )
    assert model_nsp1.gbt_compatible is True
    assert model_nsp2.soc_enabled is True
    assert model_nsp4.soc_enabled is True
    assert model_nsp4.spin_representation == "noncollinear_spinor"
