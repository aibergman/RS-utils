"""HT-15 ASE boundary and setup-source coverage."""

from __future__ import annotations

from pathlib import Path
import sys

import f90nml
import numpy as np
import pytest

from rsutil.ase import from_ase, run_setup_from_ase, to_ase, write_setup_source
from rsutil.core import OperationStatus, RSStructure, SetupResult, SetupSpec
from rsutil.setup import SetupService


ase = pytest.importorskip("ase")
from ase import Atoms  # noqa: E402


def _atoms(cell, pbc):
    atoms = Atoms(
        symbols=["Fe", "Pt", "N"],
        positions=[[0.1, 0.2, 0.3], [1.1, 1.2, 1.3], [2.1, 0.4, 1.5]],
        cell=cell,
        pbc=pbc,
    )
    atoms.set_tags([4, 7, 9])
    atoms.info["database_id"] = "example-17"
    atoms.info["rsutil_rs_labels"] = ["Fe1", "Pt1", "N1"]
    atoms.info["rsutil_rs_type_indices"] = [1, 2, 3]
    atoms.new_array("site_weights", np.asarray([1.0, 2.0, 3.0]))
    return atoms


@pytest.mark.parametrize(
    "cell",
    [
        np.diag([3.0, 4.0, 5.0]),
        np.asarray([[3.0, 0.2, 0.1], [0.0, 4.0, 0.3], [0.0, 0.0, 5.0]]),
    ],
)
@pytest.mark.parametrize("pbc", [(True, True, True), (True, True, False), (False, False, False)])
def test_ase_round_trip_preserves_geometry_order_pbc_and_metadata(cell, pbc) -> None:
    atoms = _atoms(cell, pbc)

    structure = from_ase(atoms)
    restored = to_ase(structure)

    assert structure.length_unit == "angstrom"
    assert structure.coordinate_convention.value == "cartesian"
    assert structure.pbc == pbc
    assert restored.get_chemical_symbols() == atoms.get_chemical_symbols()
    assert restored.get_atomic_numbers().tolist() == atoms.get_atomic_numbers().tolist()
    assert restored.get_pbc().tolist() == list(pbc)
    np.testing.assert_allclose(restored.cell.array, atoms.cell.array)
    np.testing.assert_allclose(restored.get_positions(), atoms.get_positions())
    np.testing.assert_array_equal(restored.get_tags(), atoms.get_tags())
    np.testing.assert_allclose(restored.arrays["site_weights"], atoms.arrays["site_weights"])
    assert restored.info["database_id"] == "example-17"
    assert restored.info["rsutil_rs_labels"] == ["Fe1", "Pt1", "N1"]
    assert restored.info["rsutil_rs_type_indices"] == [1, 2, 3]
    assert [site.rs_type_index for site in structure.sites] == [1, 2, 3]


def test_fractional_structure_is_converted_using_row_cell_orientation() -> None:
    atoms = _atoms(np.asarray([[3.0, 0.2, 0.0], [0.0, 4.0, 0.1], [0.0, 0.0, 5.0]]), (True, True, True))
    structure = from_ase(atoms)
    fractional = structure.__class__(
        lattice_vectors=structure.lattice_vectors,
        length_unit="angstrom",
        coordinate_convention="fractional",
        pbc=structure.pbc,
        sites=tuple(
            site.__class__(
                index=site.index,
                position=tuple(np.linalg.solve(np.asarray(structure.lattice_vectors).T, np.asarray(site.position))),
                element=site.element,
                atomic_number=site.atomic_number,
            )
            for site in structure.sites
        ),
    )
    restored = to_ase(fractional)
    np.testing.assert_allclose(restored.get_positions(), atoms.get_positions())
    np.testing.assert_allclose(restored.cell.array, atoms.cell.array)


def test_core_convenience_methods_use_the_optional_boundary_adapter() -> None:
    atoms = _atoms(np.diag([3.0, 4.0, 5.0]), (True, True, True))
    converted = RSStructure.from_ase(atoms)
    assert converted.to_ase().get_chemical_symbols() == atoms.get_chemical_symbols()


def test_missing_ase_has_a_clear_optional_dependency_error(monkeypatch: pytest.MonkeyPatch) -> None:
    from rsutil.ase import ASEDependencyError

    monkeypatch.setitem(sys.modules, "ase", None)
    with pytest.raises(ASEDependencyError, match="optional 'ase'.*pip install"):
        from_ase(object())


def test_setup_source_uses_existing_namelist_contract_and_preserves_species_order(tmp_path: Path) -> None:
    atoms = _atoms(np.diag([3.0, 4.0, 5.0]), (True, True, False))
    with pytest.raises(ValueError, match="experimental"):
        write_setup_source(atoms, tmp_path / "not-approved")
    bundle = write_setup_source(atoms, tmp_path / "source", experimental=True)

    lattice = f90nml.read(bundle.lattice_path)["lattice"]
    input_nml = f90nml.read(bundle.input_path)
    assert lattice["nbulk_bulk"] == 3
    assert lattice["izp"] == [1, 2, 3]
    assert np.asarray(lattice["crd"]).shape == (3, 3)
    assert input_nml["lattice"]["b3"] is False
    assert input_nml["atoms"]["label"] == ["Fe", "Pt", "N"]
    expected_wav = (3.0 * (3.0 * 4.0 * 5.0 / 3.0) / (4.0 * np.pi)) ** (1.0 / 3.0)
    np.testing.assert_allclose(input_nml["lattice"]["wav"], expected_wav)


class _ResourceAwareRunner:
    def run(self, lattice_path, input_path, spec, *, workspace, resources):
        assert tuple(Path(item).name for item in resources) == ("Fe.nml",)
        workspace = Path(workspace)
        workspace.mkdir(parents=True, exist_ok=True)
        (workspace / "lattice_new.nml").write_text("&lattice\n ntot=1\n/\n", encoding="utf-8")
        (workspace / "input_new.nml").write_text("&atoms\n label(1)='Fe'\n/\n", encoding="utf-8")
        (workspace / "radii.out").write_text("1\n1.0\n", encoding="utf-8")
        return SetupResult(
            status=OperationStatus.SUCCESS,
            candidate_files=("input_new.nml", "lattice_new.nml", "radii.out"),
        )


def test_run_setup_from_ase_delegates_to_setup_service_with_explicit_resources(tmp_path: Path) -> None:
    atoms = _atoms(np.diag([3.0, 4.0, 5.0]), (True, True, True))
    resource = tmp_path / "Fe.nml"
    resource.write_text("potential", encoding="utf-8")
    result = run_setup_from_ase(
        atoms,
        tmp_path / "destination",
        SetupSpec(),
        resources=[resource],
        experimental=True,
        setup_service=SetupService(runner=_ResourceAwareRunner()),
    )
    assert result.success
