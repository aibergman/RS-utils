"""
RS-LMTO/POSCAR/XYZ -> LMTO47 CTRL file exporter

Converts crystal structures from various formats into LMTO47 CTRL input files.
Supports RS-LMTO namelists (lattice.nml + input.nml), POSCAR, and XYZ formats.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import numpy as np

try:
    from pymatgen.core import Element
except ImportError:
    Element = None  # type: ignore

try:
    from rsutil.viz.inputreader import InputReader
    from rsutil.viz.structure import Structure
except (ImportError, AttributeError):
    InputReader = None  # type: ignore
    Structure = None  # type: ignore


ANGSTROM_TO_BOHR = 1.0 / 0.529177


def _extract_alat_from_input(input_nml_path: str) -> Optional[float]:
    try:
        import f90nml
        nml = f90nml.read(input_nml_path)
        lattice_section = nml.get("lattice", {})
        alat = lattice_section.get("alat")
        if alat is not None:
            return float(alat)
    except Exception:
        pass
    return None


def _get_atomic_number(symbol: str) -> int:
    if Element is not None:
        try:
            return Element(symbol).Z
        except Exception:
            pass

    element_z = {
        "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
        "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
        "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
        "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
        "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
        "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
        "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
        "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
        "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74,
        "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
        "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
        "Pa": 91, "U": 92, "Np": 93, "Pu": 94,
    }
    return element_z.get(symbol, 0)


def build_ctrl_payload(structure: Structure, alat_angstrom: Optional[float] = None,
                       nspin: int = 2, lmx: Optional[int] = None) -> Dict:
    lattice_angstrom = structure.unit_cell
    if lattice_angstrom is None:
        raise ValueError("Structure must have unit_cell defined")

    if alat_angstrom is None:
        alat_angstrom = np.linalg.norm(lattice_angstrom[0])

    alat_bohr = alat_angstrom * ANGSTROM_TO_BOHR
    plat = (lattice_angstrom / alat_angstrom).T

    classes = []
    species_set = set()
    for atom in structure.atoms_uc:
        if atom.species not in species_set:
            z = _get_atomic_number(atom.species)
            cls_dict = {"symbol": atom.species, "Z": z}
            if lmx is not None:
                cls_dict["lmx"] = lmx
            classes.append(cls_dict)
            species_set.add(atom.species)

    sites = []
    for atom in structure.atoms_uc:
        cart_pos_in_alat = atom.position / alat_angstrom
        sites.append({
            "symbol": atom.species,
            "pos": cart_pos_in_alat.tolist()
        })

    return {
        "alat_bohr": alat_bohr,
        "plat": plat,
        "classes": classes,
        "sites": sites,
        "nspin": nspin,
    }


def write_ctrl(payload: Dict, out_path: Union[str, Path], header: str = "RS-tools export") -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    alat_bohr = payload["alat_bohr"]
    plat = payload["plat"]
    classes = payload["classes"]
    sites = payload["sites"]

    lines = []
    lines.append(f"HEADER {header}")
    lines.append("VERS      LMASA-47")
    lines.append("IO        VERBOS=90 HELP=F WKP=F IACTIV=F ERRTOL=2 OUTPUT=LM ERR=ERR")
    lines.append(f"STRUC ALAT={alat_bohr:.8f}")

    for i, vec in enumerate(plat):
        if i == 0:
            lines.append(f"      PLAT={vec[0]:.8f} {vec[1]:.8f} {vec[2]:.8f}")
        else:
            lines.append(f"           {vec[0]:.8f} {vec[1]:.8f} {vec[2]:.8f}")

    lines.append(f"OPTIONS   NSPIN={payload['nspin']} REL=T CCOR=T")
    lines.append("SCALE     SCLWSR=T OMMAX1=.16 .18 .20 OMMAX2=.40 .45 .50")
    lines.append("START     NIT=30 BROY=T WC=-1 NMIX=1 BETA=.5")
    lines.append("          FREE=F CNVG=.0000001 CNVGET=.000001 BEGMOM=T CNTROL=T")

    for i, cls in enumerate(classes):
        cls_line = f"ATOM={cls['symbol']} Z={cls['Z']}"
        if "lmx" in cls:
            cls_line += f" LMX={cls['lmx']}"
        if i == 0:
            lines.append(f"CLASS {cls_line}")
        else:
            lines.append(f"      {cls_line}")

    for i, site in enumerate(sites):
        x, y, z = site["pos"]
        if i == 0:
            lines.append(f"SITE  ATOM={site['symbol']} POS= {x:.8f} {y:.8f} {z:.8f}")
        else:
            lines.append(f"      ATOM={site['symbol']} POS= {x:.8f} {y:.8f} {z:.8f}")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def export_ctrl(
    source: Union[Tuple[str, ...], Structure],
    out_dir: str = ".",
    header: str = "RS-tools export",
    nspin: int = 2,
    lmx: Optional[int] = None,
    **kwargs
) -> str:
    alat_angstrom = None

    if isinstance(source, tuple):
        fmt = source[0]
        if fmt == "rs_lmto":
            if len(source) < 3:
                raise ValueError("rs_lmto requires (source_type, lattice_nml_path, input_nml_path)")
            lattice_nml, input_nml = source[1], source[2]
            if InputReader is None:
                raise ImportError("InputReader not available")
            structure = InputReader.read_rs_lmto(lattice_nml, input_nml)
            alat_angstrom = _extract_alat_from_input(input_nml)
        elif fmt == "poscar":
            if len(source) < 2:
                raise ValueError("poscar requires (source_type, poscar_path)")
            poscar_path = source[1]
            if InputReader is None:
                raise ImportError("InputReader not available")
            structure = InputReader.read_poscar(poscar_path)
        elif fmt == "xyz":
            if len(source) < 2:
                raise ValueError("xyz requires (source_type, xyz_path)")
            xyz_path = source[1]
            if InputReader is None:
                raise ImportError("InputReader not available")
            structure = InputReader.read_xyz(xyz_path)
        else:
            raise ValueError(f"Unknown source format: {fmt}")
    else:
        structure = source

    if structure is None:
        raise ValueError("Failed to parse input structure")

    payload = build_ctrl_payload(structure, alat_angstrom=alat_angstrom, nspin=nspin, lmx=lmx)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ctrl_path = out_dir / "CTRL"

    write_ctrl(payload, ctrl_path, header=header)

    return str(ctrl_path)
