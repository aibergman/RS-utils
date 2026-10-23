"""Semantic regression coverage for the legacy ``rsutil`` flag workflows."""

from __future__ import annotations

import importlib
import json
import math
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path

import f90nml
import pytest

from .conftest import assert_cli_success, copy_example_files, run_rsutil


def _copy_setup_fixture(destination: Path, example: str = "bccFe") -> None:
    if example == "bccFe":
        copy_example_files(destination, example, ["Dist/lattice.nml", "Dist/input.nml"])
    else:
        copy_example_files(destination, example, ["lattice.nml", "input.nml"])


def _assert_semantically_equal(left, right, path: str = "root") -> None:
    """Compare f90nml values without making formatting/order part of the contract."""

    if isinstance(left, bool) or isinstance(right, bool):
        assert left == right, f"{path}: {left!r} != {right!r}"
    elif isinstance(left, (int, float)) and isinstance(right, (int, float)):
        assert math.isclose(float(left), float(right), rel_tol=1e-10, abs_tol=1e-10), (
            f"{path}: {left!r} != {right!r}"
        )
    elif isinstance(left, Mapping) and isinstance(right, Mapping):
        assert set(left) == set(right), f"{path}: namelist keys differ"
        for key in left:
            _assert_semantically_equal(left[key], right[key], f"{path}.{key}")
    elif isinstance(left, Sequence) and not isinstance(left, (str, bytes)) and isinstance(
        right, Sequence
    ) and not isinstance(right, (str, bytes)):
        assert len(left) == len(right), f"{path}: sequence lengths differ"
        for index, (left_item, right_item) in enumerate(zip(left, right)):
            _assert_semantically_equal(left_item, right_item, f"{path}[{index}]")
    else:
        assert left == right, f"{path}: {left!r} != {right!r}"


def _read_namelist(path: Path) -> f90nml.Namelist:
    return f90nml.read(path)


def _make_source_derived_unequal_count_fixture(destination: Path) -> None:
    """Create a worker-valid lattice with an independent source count."""

    _copy_setup_fixture(destination, "FeN")
    lattice_path = destination / "lattice.nml"
    lattice_text = lattice_path.read_text(encoding="utf-8")
    lattice_text = lattice_text.replace("nbas        = 2", "nbas        = 51")
    lattice_path.write_text(lattice_text, encoding="utf-8")


def test_setup_simple_magnetic_metal_and_reports(tmp_path: Path) -> None:
    _copy_setup_fixture(tmp_path)
    original_inputs = {
        name: (tmp_path / name).read_bytes() for name in ("lattice.nml", "input.nml")
    }

    result = run_rsutil(tmp_path, "-S", "--validate", "--summary")
    assert_cli_success(result)
    assert "All validation checks passed" in result.stdout
    assert "Summary of Changes:" in result.stdout

    lattice = _read_namelist(tmp_path / "lattice_new.nml")["lattice"]
    input_nml = _read_namelist(tmp_path / "input_new.nml")
    assert lattice["ntot"] == 1
    assert input_nml["atoms"]["label"] == ["Fe"]
    assert input_nml["control"]["recur"] == "block"
    assert input_nml["control"]["lld"] == 21
    assert input_nml["lattice"]["strux_backend"] == "strux_lib"
    assert {name: (tmp_path / name).read_bytes() for name in original_inputs} == original_inputs


def test_setup_ks_preset_and_explicit_strux_backend(tmp_path: Path) -> None:
    _copy_setup_fixture(tmp_path)
    input_path = tmp_path / "input.nml"
    input_path.write_text(
        input_path.read_text(encoding="utf-8").replace(
            "crystal_sym = 'file'",
            "crystal_sym = 'file'\nstrux_backend = 'custom_backend'",
        ),
        encoding="utf-8",
    )
    original_input = input_path.read_bytes()

    result = run_rsutil(tmp_path, "-S", "--preset", "KS", "--validate")
    assert_cli_success(result)

    generated = _read_namelist(tmp_path / "input_new.nml")
    assert generated["self"]["use_kspace"] is True
    assert generated["lattice"]["strux_backend"] == "custom_backend"
    assert input_path.read_bytes() == original_input


def test_setup_multisite_multispecies_fixture(tmp_path: Path) -> None:
    _copy_setup_fixture(tmp_path, "FeN")

    result = run_rsutil(tmp_path, "-S", "--validate")
    assert_cli_success(result)

    lattice = _read_namelist(tmp_path / "lattice_new.nml")["lattice"]
    input_nml = _read_namelist(tmp_path / "input_new.nml")
    assert lattice["ntot"] == 2
    assert input_nml["atoms"]["label"] == ["Fe1", "N1"]
    assert input_nml["lattice"]["ntype"] == 2
    assert (tmp_path / "Fe1.nml").is_file()
    assert (tmp_path / "N1.nml").is_file()


def test_setup_supercell_is_non_destructive_and_semantic(tmp_path: Path) -> None:
    _copy_setup_fixture(tmp_path)
    original_inputs = {
        name: (tmp_path / name).read_bytes() for name in ("lattice.nml", "input.nml")
    }

    result = run_rsutil(tmp_path, "-S", "--supercell", "1", "1", "2", "--validate")
    assert_cli_success(result)

    lattice = _read_namelist(tmp_path / "lattice_new.nml")["lattice"]
    input_nml = _read_namelist(tmp_path / "input_new.nml")
    assert lattice["ntot"] == 2
    assert lattice["nbas"] == 2
    assert input_nml["atoms"]["label"] == ["Fe1", "Fe2"]
    assert not (tmp_path / ".rsutil_supercell_lattice.nml").exists()
    assert not (tmp_path / ".rsutil_supercell_input.nml").exists()
    assert {name: (tmp_path / name).read_bytes() for name in original_inputs} == original_inputs


def test_setup_in_place_preserves_cli_backups(tmp_path: Path) -> None:
    _copy_setup_fixture(tmp_path)
    original_inputs = {
        name: (tmp_path / name).read_bytes() for name in ("lattice.nml", "input.nml")
    }

    result = run_rsutil(tmp_path, "-S", "--in-place", "--validate")

    assert_cli_success(result)
    assert (tmp_path / "lattice.nml.bak").read_bytes() == original_inputs["lattice.nml"]
    assert (tmp_path / "input.nml.bak").read_bytes() == original_inputs["input.nml"]
    assert not (tmp_path / "lattice_new.nml").exists()
    assert not (tmp_path / "input_new.nml").exists()


def test_seeded_sqs_setup_is_reproducible_in_isolated_runs(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    _copy_setup_fixture(first)
    _copy_setup_fixture(second)
    arguments = (
        "-S",
        "--supercell",
        "2",
        "1",
        "1",
        "--sqs",
        "Fe=0.5,Co=0.5",
        "--sqs-steps",
        "10",
        "--sqs-seed",
        "42",
    )

    assert_cli_success(run_rsutil(first, *arguments))
    assert_cli_success(run_rsutil(second, *arguments))
    _assert_semantically_equal(
        _read_namelist(first / "lattice_new.nml"),
        _read_namelist(second / "lattice_new.nml"),
    )
    _assert_semantically_equal(
        _read_namelist(first / "input_new.nml"),
        _read_namelist(second / "input_new.nml"),
    )


@pytest.mark.parametrize("example", ["bccFe", "FeN"])
def test_identical_setup_is_a_semantic_fixed_point_after_promoting_outputs(
    tmp_path: Path, example: str
) -> None:
    _copy_setup_fixture(tmp_path, example)
    arguments = ("-S", "--seed", "42")

    assert_cli_success(run_rsutil(tmp_path, *arguments))
    first_generated = {
        path.name: _read_namelist(path)
        for path in tmp_path.glob("*.nml")
        if path.name not in {"lattice.nml", "input.nml"}
    }
    assert {"lattice_new.nml", "input_new.nml"} <= set(first_generated)

    shutil.copy2(tmp_path / "lattice_new.nml", tmp_path / "lattice.nml")
    shutil.copy2(tmp_path / "input_new.nml", tmp_path / "input.nml")
    assert_cli_success(run_rsutil(tmp_path, *arguments))
    second_generated = {
        path.name: _read_namelist(path)
        for path in tmp_path.glob("*.nml")
        if path.name not in {"lattice.nml", "input.nml"}
    }

    assert set(second_generated) == set(first_generated)
    for name in first_generated:
        _assert_semantically_equal(first_generated[name], second_generated[name], path=name)


def test_legacy_setup_accepts_source_derived_unequal_counts(tmp_path: Path) -> None:
    _make_source_derived_unequal_count_fixture(tmp_path)

    source = _read_namelist(tmp_path / "lattice.nml")["lattice"]
    assert source["nbas"] == 51
    assert source["ntot"] == 2

    result = run_rsutil(tmp_path, "-S", "--validate")

    assert_cli_success(result)
    generated = _read_namelist(tmp_path / "lattice_new.nml")["lattice"]
    assert generated["ntot"] == 2
    assert generated["nbas"] == 2
    assert _read_namelist(tmp_path / "input_new.nml")["atoms"]["label"] == ["Fe1", "N1"]


def _copy_exchange_fixture(destination: Path) -> None:
    copy_example_files(
        destination,
        "FePt",
        [
            "input.nml",
            "lattice.nml",
            "clust",
            "report.out",
            "Fe.nml",
            "Pt.nml",
            "Fe_out.nml",
            "Pt_out.nml",
        ],
    )


def test_exchange_default_directories_fermi_propagation_and_rij_output(
    tmp_path: Path,
) -> None:
    _copy_exchange_fixture(tmp_path)
    original_inputs = {
        name: (tmp_path / name).read_bytes() for name in ("lattice.nml", "input.nml")
    }

    result = run_rsutil(tmp_path, "-J", "--cutoff", "1.0", "--rij-out", "rij.json")
    assert_cli_success(result)

    for atom_index in (1, 2):
        atom_dir = tmp_path / f"Atom_{atom_index}"
        assert (atom_dir / "lattice.nml").is_file()
        assert (atom_dir / "Fe_out.nml").is_file()
        assert (atom_dir / "Pt_out.nml").is_file()
        generated = _read_namelist(atom_dir / "input.nml")
        assert generated["calculation"]["post_processing"] == "exchange"
        assert math.isclose(generated["energy"]["fermi"], -0.159456, abs_tol=1e-12)
        assert generated["lattice"]["njij"] > 0
        assert len(generated["lattice"]["ijpair"]) == 2

    mappings = json.loads((tmp_path / "rij.json").read_text(encoding="utf-8"))
    assert mappings
    assert {"i", "j", "R_ij", "distance", "species_i", "species_j"} <= set(mappings[0])
    assert len(mappings[0]["R_ij"]) == 3
    assert {name: (tmp_path / name).read_bytes() for name in original_inputs} == original_inputs


def test_exchange_unique_soc_interface_is_non_crashing(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)

    result = run_rsutil(tmp_path, "-J", "--cutoff", "1.0", "--unique", "--soc")
    assert_cli_success(result)
    assert (tmp_path / "Atom_1" / "input.nml").is_file()
    assert (tmp_path / "Atom_2" / "input.nml").is_file()


def test_legacy_exchange_explicit_supported_reduction_is_semantic(tmp_path: Path) -> None:
    full_root = tmp_path / "full"
    reduced_root = tmp_path / "reduced"
    full_root.mkdir()
    reduced_root.mkdir()
    _copy_exchange_fixture(full_root)
    _copy_exchange_fixture(reduced_root)

    full_result = run_rsutil(full_root, "-J", "--cutoff", "1.0", "--full-pairs")
    reduced_result = run_rsutil(
        reduced_root,
        "-J",
        "--cutoff",
        "1.0",
        "--unique",
        "--soc",
        "--magnetic-reference",
        "collinear_fm",
        "--magnetic-axis",
        "0",
        "0",
        "1",
    )

    assert_cli_success(full_result)
    assert_cli_success(reduced_result)
    full_jobs = sorted(full_root.glob("Atom_*/input.nml"))
    reduced_jobs = sorted(reduced_root.glob("Atom_*/input.nml"))
    assert reduced_jobs
    full_pair_count = sum(int(_read_namelist(path)["lattice"]["njij"]) for path in full_jobs)
    reduced_pair_count = sum(int(_read_namelist(path)["lattice"]["njij"]) for path in reduced_jobs)
    assert reduced_pair_count < full_pair_count
    for path in reduced_jobs:
        generated = _read_namelist(path)
        assert generated["calculation"]["post_processing"] == "exchange"
        assert generated["lattice"]["njij"] > 0
        manifest = json.loads((path.parent / "rsutil_exchange_manifest.json").read_text(encoding="utf-8"))
        assert manifest["spec"]["magnetic_reference"] == "collinear_fm"
        assert manifest["spec"]["symmetry_reduction_applied"] is True


def test_existing_band_and_export_modules_import() -> None:
    for module_name in (
        "rsutil.bands",
        "rsutil.bands.readers.rs_lmto",
        "rsutil.export",
        "rsutil.export.lmto47",
        "rsutil.export.uppasd",
    ):
        assert importlib.import_module(module_name) is not None
