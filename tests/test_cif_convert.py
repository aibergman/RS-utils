"""CIF -> RS-LMTO source namelist conversion, checked against cif2cell output."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

f90nml = pytest.importorskip("f90nml")
pytest.importorskip("pymatgen")

from rsutil.cif import CifConversionError, convert_cif, read_cif, write_rs_input

DATA = Path(__file__).parent / "data" / "cif"
CIF = DATA / "Ni2MnGa_premartensic.cif"
REF_LATTICE = DATA / "Ni2MnGa_premartensic.cif2cell.lattice.nml"
REF_INPUT = DATA / "Ni2MnGa_premartensic.cif2cell.input.nml"

BCC_FE = """data_Fe
_cell_length_a 2.87
_cell_length_b 2.87
_cell_length_c 2.87
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'I m -3 m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0 0 0 1.0
"""

DISORDERED = BCC_FE.replace(
    "_atom_site_occupancy\nFe1 Fe 0 0 0 1.0\n",
    "_atom_site_occupancy\nFe1 Fe 0 0 0 0.5\nCo1 Co 0 0 0 0.5\n",
)


def _crd_by_element(lattice, inp):
    crd = np.asarray(lattice["crd"], dtype=float)
    crd = crd if crd.shape[1] == 3 else crd.T
    grouped: dict[str, list] = {}
    for label, row in zip(inp["atoms"]["label"], crd):
        grouped.setdefault(label.rstrip("0123456789"), []).append(tuple(np.round(row, 6)))
    return {key: np.array(sorted(rows)) for key, rows in grouped.items()}


@pytest.fixture(scope="module")
def converted(tmp_path_factory):
    out = tmp_path_factory.mktemp("ni2mnga")
    files, conversion = convert_cif(CIF, out, ct=4.30)
    return files, conversion


def test_matches_cif2cell_reference(converted):
    files, _ = converted
    new_l = f90nml.read(str(files.lattice))["lattice"]
    ref_l = f90nml.read(str(REF_LATTICE))["lattice"]
    new_i = f90nml.read(str(files.input))
    ref_i = f90nml.read(str(REF_INPUT))

    for key in ("nbulk_bulk", "ntot", "nbas", "nrec"):
        assert new_l[key] == ref_l[key] == 24
    for key in ("izp", "no", "iu", "ib", "irec"):
        assert list(new_l[key]) == list(ref_l[key]) == list(range(1, 25))
    assert np.allclose(new_l["a"], ref_l["a"], atol=1e-9)

    # Atom order within an element is a cif2cell implementation detail, so
    # positions are compared per element; labels and ordering of elements must match.
    assert list(new_i["atoms"]["label"]) == list(ref_i["atoms"]["label"])
    new_crd, ref_crd = _crd_by_element(new_l, new_i), _crd_by_element(ref_l, ref_i)
    assert new_crd.keys() == ref_crd.keys()
    for element in ref_crd:
        assert np.allclose(new_crd[element], ref_crd[element], atol=1e-6)

    for key in ("rc", "alat", "crystal_sym", "ntype", "r2"):
        assert new_i["lattice"][key] == ref_i["lattice"][key]
    assert list(new_i["lattice"]["ct"]) == list(ref_i["lattice"]["ct"])
    for group in ("calculation", "self", "energy", "control", "mix", "hamiltonian"):
        assert dict(new_i[group]) == dict(ref_i[group])


def test_structure_contract(converted):
    _, conversion = converted
    structure = conversion.structure
    assert structure.length_unit == "angstrom"
    assert conversion.alat == pytest.approx(4.116)
    assert [site.rs_type_index for site in structure.sites] == list(range(1, 25))
    assert [site.element for site in structure.sites] == ["Ga"] * 6 + ["Ni"] * 12 + ["Mn"] * 6


def test_template_group_order_preserved(converted):
    files, _ = converted
    text = files.input.read_text()
    groups = [line for line in text.splitlines() if line.startswith("&")]
    assert groups == [
        "&lattice", "&atoms", "&calculation", "&self",
        "&energy", "&control", "&mix", "&hamiltonian",
    ]


def test_template_and_wav_override(tmp_path):
    template = tmp_path / "template.nml"
    template.write_text("&control\ncalctype='B'\nnsp=2\n/\n&lattice\nwav=9.9\nrc=7\n/\n")
    out = tmp_path / "out"
    files, _ = convert_cif(
        _write(tmp_path, BCC_FE), out, ct=3.3, template=template, wav=1.4088
    )
    nml = f90nml.read(str(files.input))
    assert nml["control"]["nsp"] == 2
    assert "calculation" not in nml  # only the template's groups are used
    assert nml["lattice"]["wav"] == pytest.approx(1.4088)  # explicit wins over template
    assert nml["lattice"]["rc"] == 50  # generated value wins over template


def _write(tmp_path: Path, text: str, name: str = "x.cif") -> Path:
    path = tmp_path / name
    path.write_text(text)
    return path


def test_bcc_unit_cell(tmp_path):
    files, conversion = convert_cif(_write(tmp_path, BCC_FE), tmp_path / "o", ct=3.3)
    lattice = f90nml.read(str(files.lattice))["lattice"]
    assert lattice["ntot"] == 2
    assert np.allclose(np.asarray(lattice["a"]).reshape(3, 3), np.eye(3))
    crd = np.asarray(lattice["crd"], dtype=float)
    crd = crd if crd.shape[1] == 3 else crd.T
    assert np.allclose(sorted(map(tuple, crd)), [(0, 0, 0), (0.5, 0.5, 0.5)])
    assert conversion.alat == pytest.approx(2.87)


def test_supercell_keeps_primitive_alat(tmp_path):
    files, conversion = convert_cif(
        _write(tmp_path, BCC_FE), tmp_path / "o", ct=3.3, supercell=(2, 1, 1)
    )
    lattice = f90nml.read(str(files.lattice))["lattice"]
    assert lattice["ntot"] == 4
    assert conversion.alat == pytest.approx(2.87)
    # f90nml returns a(:, i) as the i-th column of a 3x3 array: a1 = (2, 0, 0) alat.
    a = np.asarray(lattice["a"], dtype=float)
    assert a.shape == (3, 3)
    assert np.allclose(a, np.diag([2.0, 1.0, 1.0]))


def test_rejects_partial_occupancy(tmp_path):
    with pytest.raises(CifConversionError, match="partially occupied"):
        read_cif(_write(tmp_path, DISORDERED))


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
def test_rejects_bad_ct(tmp_path, bad):
    conversion = read_cif(_write(tmp_path, BCC_FE))
    out = tmp_path / "o"
    with pytest.raises(CifConversionError, match="ct"):
        write_rs_input(conversion, out, ct=bad)
    assert not out.exists()  # nothing written on validation failure


def test_rejects_bad_supercell_and_alat_and_missing_file(tmp_path):
    cif = _write(tmp_path, BCC_FE)
    with pytest.raises(CifConversionError, match="supercell"):
        read_cif(cif, supercell=(0, 1, 1))
    with pytest.raises(CifConversionError, match="alat"):
        read_cif(cif, alat=-1)
    with pytest.raises(CifConversionError, match="does not exist"):
        read_cif(tmp_path / "missing.cif")
    with pytest.raises(CifConversionError):
        read_cif(_write(tmp_path, "not a cif at all", "bad.cif"))


def test_refuses_to_overwrite(tmp_path):
    cif = _write(tmp_path, BCC_FE)
    out = tmp_path / "o"
    convert_cif(cif, out, ct=3.3)
    with pytest.raises(FileExistsError):
        convert_cif(cif, out, ct=3.3)
    convert_cif(cif, out, ct=3.3, overwrite=True)


def test_cli_setup_cif(tmp_path):
    out = tmp_path / "run"
    cmd = [sys.executable, "-m", "rsutil.cli", "setup", "cif", str(CIF), "--ct", "4.3", "-o", str(out)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert (out / "lattice.nml").is_file() and (out / "input.nml").is_file()
    assert "Ga6, Ni12, Mn6" in result.stdout

    again = subprocess.run(cmd, capture_output=True, text=True)
    assert again.returncode == 1 and "overwrite" in again.stderr
