"""CIF -> RS-LMTO source namelist conversion, checked against cif2cell output."""

from __future__ import annotations

import os
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


def test_bcc_conventional_cell(tmp_path):
    files, conversion = convert_cif(
        _write(tmp_path, BCC_FE), tmp_path / "o", ct=3.3, cell="conventional"
    )
    lattice = f90nml.read(str(files.lattice))["lattice"]
    assert lattice["ntot"] == 2
    assert np.allclose(np.asarray(lattice["a"]).reshape(3, 3), np.eye(3))
    crd = np.asarray(lattice["crd"], dtype=float)
    crd = crd if crd.shape[1] == 3 else crd.T
    assert np.allclose(sorted(map(tuple, crd)), [(0, 0, 0), (0.5, 0.5, 0.5)])
    assert conversion.alat == pytest.approx(2.87)


def test_supercell_keeps_conventional_alat(tmp_path):
    files, conversion = convert_cif(
        _write(tmp_path, BCC_FE),
        tmp_path / "o",
        ct=3.3,
        supercell=(2, 1, 1),
        cell="conventional",
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


# --- automatic ct ----------------------------------------------------------

from rsutil.cif import neighbor_counts, suggest_ct
from rsutil.core import RSStructure, StructureSite


def _structure(cell, positions, element="Fe"):
    return RSStructure(
        lattice_vectors=np.asarray(cell, dtype=float),
        length_unit="angstrom",
        sites=tuple(
            StructureSite(index=i, position=tuple(p), element=element)
            for i, p in enumerate(positions)
        ),
    )


def test_suggest_ct_bcc_gives_15_atoms():
    a = 2.87
    bcc = _structure(np.eye(3) * a, [(0, 0, 0), (a / 2, a / 2, a / 2)])
    ct = suggest_ct(bcc)
    assert set(neighbor_counts(bcc, ct)) == {15}
    assert a < ct < a * np.sqrt(2)  # between the 2nd and 3rd shells


def test_suggest_ct_fcc_gives_19_atoms():
    a = 3.6
    fcc = _structure(
        np.eye(3) * a, [(0, 0, 0), (a / 2, a / 2, 0), (a / 2, 0, a / 2), (0, a / 2, a / 2)], "Cu"
    )
    assert set(neighbor_counts(fcc, suggest_ct(fcc))) == {19}


def test_suggest_ct_reproduces_feco_reference_cutoff():
    feco = _structure(np.eye(3) * 3.0, [(0, 0, 0), (1.5, 1.5, 1.5)])
    assert suggest_ct(feco) == pytest.approx(3.3)


def test_suggest_ct_respects_min_atoms_and_shell_splitting(converted):
    _, conversion = converted
    structure = conversion.structure  # distorted, shells split by ~0.05 A
    for wanted in (15, 19, 27):
        ct = suggest_ct(structure, min_atoms=wanted)
        counts = neighbor_counts(structure, ct)
        assert min(counts) >= wanted
        # no site may gain an atom from a split shell: counts are uniform
        assert len(set(counts)) == 1


@pytest.mark.parametrize("bad", [1, 0, -3, True])
def test_suggest_ct_rejects_bad_min_atoms(bad):
    with pytest.raises(CifConversionError, match="ct_neighbors"):
        suggest_ct(_structure(np.eye(3) * 2.87, [(0, 0, 0)]), min_atoms=bad)


def test_ct_defaults_automatically_and_is_reported(tmp_path):
    files, _ = convert_cif(_write(tmp_path, BCC_FE), tmp_path / "o", cell="conventional")
    assert files.ct_auto and files.neighbors_min == files.neighbors_max == 15
    nml = f90nml.read(str(files.input))["lattice"]
    assert list(nml["ct"]) == [pytest.approx(files.ct, abs=0.005)] * 2
    assert nml["r2"] == pytest.approx(files.ct**2)
    explicit, _ = convert_cif(_write(tmp_path, BCC_FE), tmp_path / "p", ct=3.3)
    assert not explicit.ct_auto and explicit.ct == 3.3


# --- setup bulk --cif ------------------------------------------------------


def _run_bulk(cwd, *extra):
    cmd = [sys.executable, "-m", "rsutil.cli", "setup", "bulk", *extra]
    # setup bulk works in the current directory, so run there and make the
    # (possibly uninstalled) package importable from the repository root.
    env = {**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1])}
    return subprocess.run(cmd, capture_output=True, text=True, cwd=cwd, env=env)


def test_setup_bulk_cif_writes_sources_before_setup(tmp_path):
    cif = _write(tmp_path, BCC_FE, "fe.cif")
    work = tmp_path / "work"
    work.mkdir()
    result = _run_bulk(work, "--cif", str(cif))
    # Whether the Fortran kernel is available decides the exit status; the
    # conversion step must have happened and been reported either way.
    assert "Converted" in result.stdout and "chosen automatically" in result.stdout
    assert (work / "lattice.nml").is_file() and (work / "input.nml").is_file()
    # body-centred Fe: the default is the 1-atom primitive cell
    assert f90nml.read(str(work / "lattice.nml"))["lattice"]["ntot"] == 1
    conventional = tmp_path / "conv"
    conventional.mkdir()
    result = _run_bulk(conventional, "--cif", str(cif), "--cell", "conventional")
    assert f90nml.read(str(conventional / "lattice.nml"))["lattice"]["ntot"] == 2


def test_setup_bulk_cif_refuses_to_overwrite_sources(tmp_path):
    cif = _write(tmp_path, BCC_FE, "fe.cif")
    work = tmp_path / "work"
    work.mkdir()
    (work / "lattice.nml").write_text("&lattice\n/\n")
    result = _run_bulk(work, "--cif", str(cif))
    assert result.returncode == 1 and "--overwrite-source" in result.stderr
    assert (work / "lattice.nml").read_text() == "&lattice\n/\n"  # untouched


def test_setup_bulk_cif_rejects_disorder(tmp_path):
    cif = _write(tmp_path, DISORDERED, "alloy.cif")
    work = tmp_path / "work"
    work.mkdir()
    result = _run_bulk(work, "--cif", str(cif))
    assert result.returncode == 1 and "partially occupied" in result.stderr
    assert not (work / "lattice.nml").exists()


def test_setup_bulk_without_cif_is_unchanged(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    result = _run_bulk(work)  # no sources, no --cif: old missing-input error
    assert result.returncode == 1
    assert "Converted" not in result.stdout


# --- primitive vs conventional cell ----------------------------------------

FEN = DATA / "FeN_zincblende.cif"


def _a_crd(path):
    lattice = f90nml.read(str(path))["lattice"]
    a = np.asarray(lattice["a"], dtype=float)
    crd = np.asarray(lattice["crd"], dtype=float)
    return lattice, a, (crd if crd.shape[1] == 3 else crd.T)


def test_fcc_primitive_is_default_and_matches_cif2cell(tmp_path):
    # Zinc-blende FeN: the CIF holds the 8-atom cubic cell; cif2cell gives the
    # 2-atom primitive cell with these vectors and positions (units of alat).
    files, conversion = convert_cif(FEN, tmp_path / "o", ct=3.0)
    lattice, a, crd = _a_crd(files.lattice)
    assert lattice["ntot"] == lattice["nbulk_bulk"] == lattice["nbas"] == lattice["nrec"] == 2
    assert conversion.alat == pytest.approx(4.307)  # conventional length, not |a1|
    assert np.allclose(a, [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
    assert np.allclose(crd, [[0, 0, 0], [0.25, 0.25, 0.25]])
    assert list(lattice["izp"]) == [1, 2]
    assert conversion.structure.metadata["centering"] == "F"
    assert conversion.structure.metadata["cell"] == "primitive"
    assert [site.element for site in conversion.structure.sites] == ["Fe", "N"]


def test_conventional_cell_option_keeps_cubic_cell(tmp_path):
    files, conversion = convert_cif(FEN, tmp_path / "o", ct=3.0, cell="conventional")
    lattice, a, crd = _a_crd(files.lattice)
    assert lattice["ntot"] == 8
    assert np.allclose(a, np.eye(3))
    assert conversion.structure.metadata["cell"] == "conventional"
    assert [site.element for site in conversion.structure.sites] == ["Fe"] * 4 + ["N"] * 4


def test_primitive_and_conventional_describe_the_same_crystal(tmp_path):
    """Same density and same neighbour environment, whichever cell is written."""
    prim = read_cif(FEN, cell="primitive").structure
    conv = read_cif(FEN, cell="conventional").structure
    volume = lambda s: abs(np.linalg.det(np.asarray(s.lattice_vectors)))  # noqa: E731
    assert volume(conv) / volume(prim) == pytest.approx(4.0)
    for ct in (3.0, 4.5):
        assert sorted(set(neighbor_counts(prim, ct))) == sorted(set(neighbor_counts(conv, ct)))


def test_primitive_default_leaves_primitive_lattices_untouched(converted):
    # Ni2MnGa is a P lattice: the default must not alter or rotate it.
    _, conversion = converted
    assert conversion.structure.metadata["centering"] == "P"
    assert len(conversion.structure.sites) == 24
    assert np.allclose(
        np.asarray(conversion.structure.lattice_vectors) / conversion.alat,
        np.diag([1.0, 12.35 / 4.116, 5.822 / 4.116]),
    )


CENTERED = """data_x
_cell_length_a {a}
_cell_length_b {b}
_cell_length_c {c}
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma {gamma}
_symmetry_space_group_name_H-M '{sg}'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0 0 0 1.0
"""


@pytest.mark.parametrize(
    "sg, a, b, c, gamma, letter, conv_atoms",
    [
        ("C 1 2 1", 3.0, 4.0, 5.0, 90, "C", 2),   # base-centred
        ("I 4/m m m", 3.0, 3.0, 5.0, 90, "I", 2),  # body-centred tetragonal
        ("F m m m", 3.0, 4.0, 5.0, 90, "F", 4),   # face-centred orthorhombic
        ("P m m m", 3.0, 4.0, 5.0, 90, "P", 1),   # primitive: unchanged
    ],
)
def test_centering_types_reduce_to_one_atom(tmp_path, sg, a, b, c, gamma, letter, conv_atoms):
    cif = _write(tmp_path, CENTERED.format(sg=sg, a=a, b=b, c=c, gamma=gamma), "c.cif")
    conventional = read_cif(cif, cell="conventional")
    primitive = read_cif(cif, cell="primitive")
    assert len(conventional.structure.sites) == conv_atoms
    assert len(primitive.structure.sites) == 1
    assert primitive.structure.metadata["centering"] == letter
    assert not primitive.warnings or all("primitive" not in w for w in primitive.warnings)
    volume = lambda s: abs(np.linalg.det(np.asarray(s.lattice_vectors)))  # noqa: E731
    assert volume(conventional.structure) / volume(primitive.structure) == pytest.approx(conv_atoms)
    assert primitive.alat == pytest.approx(conventional.alat)


def test_rejects_unknown_cell_choice(tmp_path):
    with pytest.raises(CifConversionError, match="cell"):
        read_cif(_write(tmp_path, BCC_FE), cell="niggli")


def test_cli_cell_option(tmp_path):
    out_p, out_c = tmp_path / "p", tmp_path / "c"
    base = [sys.executable, "-m", "rsutil.cli", "setup", "cif", str(FEN), "--ct", "3.0"]
    assert subprocess.run(base + ["-o", str(out_p)], capture_output=True, text=True).returncode == 0
    assert subprocess.run(
        base + ["-o", str(out_c), "--cell", "conventional"], capture_output=True, text=True
    ).returncode == 0
    assert f90nml.read(str(out_p / "lattice.nml"))["lattice"]["ntot"] == 2
    assert f90nml.read(str(out_c / "lattice.nml"))["lattice"]["ntot"] == 8


RHOMBOHEDRAL = """data_r
_cell_length_a 3.0
_cell_length_b 3.0
_cell_length_c 8.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
loop_
_symmetry_equiv_pos_as_xyz
'x, y, z'
'x+2/3, y+1/3, z+1/3'
'x+1/3, y+2/3, z+2/3'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0 0 0 1.0
"""


def test_rhombohedral_centering_reduces_to_one_atom(tmp_path):
    cif = _write(tmp_path, RHOMBOHEDRAL, "r.cif")
    conventional = read_cif(cif, cell="conventional")
    primitive = read_cif(cif, cell="primitive")
    assert len(conventional.structure.sites) == 3
    assert len(primitive.structure.sites) == 1
    assert primitive.structure.metadata["centering"] == "R"
    volume = lambda s: abs(np.linalg.det(np.asarray(s.lattice_vectors)))  # noqa: E731
    assert volume(conventional.structure) / volume(primitive.structure) == pytest.approx(3.0)
