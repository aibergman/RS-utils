"""
RS-LMTO -> UppASD exporter (starter module)

This module provides a minimal, extensible exporter that converts structural
and (optional) magnetic data from RS-LMTO inputs (lattice.nml, input.nml,
clust, etc.) into a text output that can be adapted for UppASD.

Important: UppASD has multiple input conventions. This file provides a
convert-and-template approach — it will write a small, well-documented
structure file (positions, species, lattice) plus an optional moments file
and a placeholder for exchange interactions. To finish the integration I need
an example of the exact UppASD input format you want to target (or a sample
UppASD input file) and some RS-LMTO example outputs to verify parsing and
unit conversions.

What I expect from you to proceed:
- Example RS-LMTO files: `lattice.nml`, `input.nml`, and optionally `clust`,
  and any generated `*_out.nml` or `report.out` files containing Fermi energy
  or magnetic moments.
- A sample UppASD structure file (or documentation) showing exact layout for
  lattice, atoms, spin magnitudes, and exchange interaction files (if you
  want Jij exported).
- Rules for mapping species -> moment magnitudes if moments are not available
  (defaults, units: Bohr magneton vs muB, lattice units).

This starter module provides:
- parse_rs_lmto_nml(lattice_nml, input_nml): read RS-LMTO via existing f90nml
  reader (relies on the project's `rsutil` Python utilities if available).
- parse_rs_clust(clust_file, input_file): convenience wrapper to read `clust`.
- build_uppasd_payload(structure, moments=None, j_matrix=None): prepare the
  data dict to be rendered.
- write_simple_uppasd_structure(payload, out_prefix, template=None): writes
  three text files that you can adapt:
    - {out_prefix}.struct (structure: lattice + atoms)
    - {out_prefix}.mom  (optional moments, one per atom)
    - {out_prefix}.jij  (optional exchange pair list placeholder)

Extend this file once you provide a target UppASD input example and sample
RS-LMTO outputs. The module intentionally keeps the exporter generic and
template-driven so we can adapt to whichever UppASD variant you use.

"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

try:
    # Prefer using the project's existing InputReader & Structure types
    from rsutil.viz.inputreader import InputReader
    from rsutil.viz.structure import Structure
except (ImportError, AttributeError):
    InputReader = None  # type: ignore
    Structure = None  # type: ignore


def parse_rs_lmto_nml(lattice_nml: str, input_nml: str) -> object:
    """
    Parse RS-LMTO files and return a Structure-like object.

    This function attempts to use the project's `InputReader.read_rs_lmto_nml`
    if available; otherwise it returns a minimal dict with lattice and atoms.

    Returns
    -------
    structure
        Either a `Structure` instance (from rsutil.viz.structure) or a dict with
        keys: unit_cell (3x3 array-like), atoms: list of (species, [x,y,z])
    """
    if InputReader is not None and hasattr(InputReader, "read_rs_lmto_nml"):
        struct = InputReader.read_rs_lmto_nml(lattice_nml, input_nml)
        return struct

    # Fallback: try parse with f90nml (best-effort)
    import f90nml

    lat = f90nml.read(lattice_nml).get("lattice", {})
    inp = f90nml.read(input_nml)

    # Attempt to extract lattice vectors and basis coords similar to InputReader
    a1 = lat.get("a(:, 1)")
    a2 = lat.get("a(:, 2)")
    a3 = lat.get("a(:, 3)")
    if a1 and a2 and a3:
        unit_cell = [a1, a2, a3]
    else:
        unit_cell = None

    crd1 = lat.get("crd(:, 1)") or []
    crd2 = lat.get("crd(:, 2)") or []
    crd3 = lat.get("crd(:, 3)") or []
    atoms = []
    n = len(crd1)
    species = []
    atoms_block = inp.get("atoms", {})
    i = 1
    # read label(i) keys
    while True:
        key = f"label({i})"
        if key in atoms_block:
            raw = atoms_block[key]
            species.append("".join(ch for ch in str(raw).strip("'\"") if ch.isalpha()))
            i += 1
        else:
            break

    if species and n > 0:
        if len(species) == 1:
            species = species * n
        for idx in range(n):
            atoms.append((species[idx], [crd1[idx], crd2[idx], crd3[idx]]))

    return {"unit_cell": unit_cell, "atoms": atoms}


def parse_rs_clust(clust_file: str, input_file: str) -> object:
    """Parse `clust` + `input.nml` pair and return structure-like object.

    Uses the project's `InputReader.read_rs_clust` when available.
    """
    if InputReader is not None and hasattr(InputReader, "read_rs_clust"):
        return InputReader.read_rs_clust(clust_file, input_file)

    # Fallback: very small parser reading species from input.nml and coords from clust
    atoms = []
    with open(clust_file, "r", encoding="utf-8") as f:
        lines = [ln.rstrip() for ln in f if ln.strip()]
    # naive parsing: skip header and read floats in groups
    if not lines:
        raise ValueError("Empty clust file")
    # get species lookup from input.nml
    import f90nml

    inp = f90nml.read(input_file)
    atoms_block = inp.get("atoms", {})
    species_lookup = []
    i = 1
    while True:
        key = f"label({i})"
        if key in atoms_block:
            raw = atoms_block[key]
            species_lookup.append("".join(ch for ch in str(raw).strip("'\"") if ch.isalpha()))
            i += 1
        else:
            break

    # Simplified: read coordinates assuming each data line contains x y z type ...
    for line in lines[1:]:
        parts = line.split()
        if len(parts) >= 4:
            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            type_idx = int(parts[3]) - 1
            species = species_lookup[type_idx] if type_idx < len(species_lookup) else "X"
            atoms.append((species, [x, y, z]))
    return {"unit_cell": None, "atoms": atoms}


def build_uppasd_payload(structure: object, moments: Optional[List[float]] = None,
                         j_pairs: Optional[List[Tuple[int, int, float]]] = None) -> Dict:
    """
    Convert a structure-like object into a payload ready to be rendered as
    UppASD inputs.

    Parameters
    ----------
    structure : Structure or dict
        If a `Structure` instance from `rsutil.viz.structure`, the exporter will
        try to use its attributes. Otherwise, expect a dict with keys:
        - unit_cell: 3x3 array-like (each row a lattice vector), or None
        - atoms: list of (species, [x,y,z])
    moments : optional list of floats
        Magnetic moment per atom (same ordering as atoms). If None, zeros
        are used or extracted if available on the structure object.
    j_pairs : optional list of (i, j, Jij)
        Pairwise exchange interactions (1-based or 0-based indices — you must
        confirm which you prefer). This will be written to a simple CSV.

    Returns
    -------
    payload : dict
    """
    payload = {"lattice": None, "atoms": [], "moments": [], "jij": []}

    # Structure is a project Structure? try duck-typing
    if Structure is not None and isinstance(structure, Structure):
        # Try to extract unit cell and atom list from the Structure object
        payload["lattice"] = getattr(structure, "unit_cell", None)
        # pick unit cell atoms (unit cell representation)
        atoms_obj = getattr(structure, "atoms_uc", None)
        if atoms_obj is None:
            # Try generic attribute
            atoms_obj = getattr(structure, "atoms", None)
        if atoms_obj:
            for atom in atoms_obj:
                # Atom objects in glyphs have .species and .position
                spec = getattr(atom, "species", getattr(atom, "atom_type", "X"))
                pos = getattr(atom, "position", None)
                payload["atoms"].append((spec, list(pos)))
    else:
        # dict fallback
        payload["lattice"] = structure.get("unit_cell") if isinstance(structure, dict) else None
        payload["atoms"] = structure.get("atoms", []) if isinstance(structure, dict) else []

    # moments
    if moments is not None:
        payload["moments"] = list(moments)
    else:
        # try to find spins on structure
        if Structure is not None and isinstance(structure, Structure):
            spins = getattr(structure, "spins_uc", [])
            if spins:
                payload["moments"] = [getattr(s, "magnitude", 0.0) for s in spins]
            else:
                payload["moments"] = [0.0] * len(payload["atoms"])
        else:
            payload["moments"] = [0.0] * len(payload["atoms"])

    if j_pairs:
        payload["jij"] = list(j_pairs)

    return payload


def write_simple_uppasd_structure(payload: Dict, out_prefix: str,
                                  template: Optional[str] = None) -> List[str]:
    """
    Write a set of simple files representing the structure in a format you can
    adapt for UppASD. Files written (with `out_prefix`):
      - {out_prefix}.struct : a small text file with lattice and atom list
      - {out_prefix}.mom    : optional moments file (one value per atom)
      - {out_prefix}.jij    : optional CSV for exchange pairs: i,j,Jij

    The actual UppASD format may differ; once you provide an example I'll
    adapt the exporter to match the canonical layout.

    Returns
    -------
    list[str]
        Paths of written files.
    """
    out_prefix = str(out_prefix)
    written = []

    struct_path = out_prefix + ".struct"
    with open(struct_path, "w", encoding="utf-8") as fh:
        fh.write("# RS-LMTO -> UppASD auto-generated structure file\n")
        fh.write("# Edit or provide a UppASD template to match your target format\n\n")
        lattice = payload.get("lattice")
        if lattice is not None:
            fh.write("# Lattice vectors (matrix rows are vectors)\n")
            for row in lattice:
                fh.write(" ".join(str(x) for x in row) + "\n")
            fh.write("\n")
        fh.write("# Atoms: species x y z\n")
        for spec, pos in payload.get("atoms", []):
            fh.write(f"{spec} {pos[0]: .12f} {pos[1]: .12f} {pos[2]: .12f}\n")
    written.append(struct_path)

    # moments
    mom_path = out_prefix + ".mom"
    moments = payload.get("moments")
    if moments:
        with open(mom_path, "w", encoding="utf-8") as fh:
            fh.write("# Magnetic moments (one per atom, same ordering as .struct)\n")
            for m in moments:
                fh.write(f"{m: .8f}\n")
        written.append(mom_path)

    # jij
    jij_path = out_prefix + ".jij"
    jlist = payload.get("jij")
    if jlist:
        with open(jij_path, "w", encoding="utf-8") as fh:
            fh.write("# i j Jij\n")
            for i, j, Jij in jlist:
                fh.write(f"{i} {j} {Jij: .8f}\n")
        written.append(jij_path)

    return written


def _write_posfile(payload: Dict, path: str) -> str:
    """Write posfile in UppASD expected format (index, type, x y z).

    Coordinates are written in Cartesian (alat=1) as requested.
    """
    atoms = payload.get("atoms", [])
    # determine type ids: map species -> type id in order of appearance
    species_to_type = {}
    types = []
    for spec, _ in atoms:
        if spec not in species_to_type:
            species_to_type[spec] = len(species_to_type) + 1
        types.append(species_to_type[spec])

    with open(path, "w", encoding="utf-8") as fh:
        for idx, ((spec, pos), t) in enumerate(zip(atoms, types), start=1):
            fh.write(f"{idx} {t}   {pos[0]: .6f}  {pos[1]: .6f}  {pos[2]: .6f}\n")
    return path


def _write_momfile(payload: Dict, path: str) -> str:
    """Write momfile: index, chemistry(=1 for ordered), magnitude (muB), ux uy uz

    Default direction is (1,0,0) if no direction is provided.
    """
    atoms = payload.get("atoms", [])
    moments = payload.get("moments", [0.0] * len(atoms))
    with open(path, "w", encoding="utf-8") as fh:
        for idx, m in enumerate(moments, start=1):
            # chemistry field: 1 for ordered structures
            chem = 1
            # direction default
            ux, uy, uz = 1.0, 0.0, 0.0
            fh.write(f"{idx} {chem} {m: .6f} {ux: .6f} {uy: .6f} {uz: .6f}\n")
    return path


def _write_jfile(payload: Dict, path: str) -> str:
    """Write jfile: type_i type_j rx ry rz Jij (Jij in mRy)

    We assume payload['jij'] is a list of (i, j, Jij) tuples where i,j are
    type ids (1-based). If rij is provided in the tuple, accept both formats.
    """
    jlist = payload.get("jij", [])
    with open(path, "w", encoding="utf-8") as fh:
        for entry in jlist:
            if len(entry) == 3:
                i, j, Jij = entry
                # default rij unknown -> write zeros
                rx, ry, rz = 0.0, 0.0, 0.0
            elif len(entry) >= 6:
                i, j, rx, ry, rz, Jij = entry[:6]
            else:
                # skip malformed
                continue
            fh.write(f"{int(i)} {int(j)}  {rx: .6f}  {ry: .6f}  {rz: .6f}  {Jij: .14f}\n")
    return path


def write_inpsd_dat(out_dir: str, simid: str = "rs_export", ncell: Tuple[int, int, int] = (1, 1, 1),
                    BC: Tuple[str, str, str] = ("P", "P", "P"), cell: Optional[object] = None,
                    Sym: int = 1, posfile: str = "./posfile", momfile: str = "./momfile",
                    exchange: str = "./jfile", dm: Optional[str] = None) -> str:
    """Write the main `inpsd.dat` file referencing the external files.

    This follows the example structure in UppASD_format.md.
    """
    path = os.path.join(out_dir, "inpsd.dat")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"simid {simid}\n")
        fh.write(f"ncell     {ncell[0]:d}        {ncell[1]:d}       {ncell[2]:d}\n")
        fh.write(f"BC        {BC[0]}         {BC[1]}         {BC[2]}\n")
        # write cell vectors; if not provided, write identity
        if cell is None:
            cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        fh.write("cell      {:.5f}   {:.5f}   {:.5f}\n".format(*cell[0]))
        fh.write("          {:.5f}   {:.5f}   {:.5f}\n".format(*cell[1]))
        fh.write("          {:.5f}   {:.5f}   {:.5f}\n".format(*cell[2]))
        fh.write(f"Sym       {Sym}\n\n")
        fh.write(f"posfile    {posfile}\n")
        fh.write(f"momfile    {momfile}\n")
        fh.write(f"exchange   {exchange}\n")
        if dm:
            fh.write(f"dm         {dm}\n")
    return path


def parse_jijso(jij_path: str):
    """Parse a `jijso.out`-style file into a list of (i,j,rx,ry,rz,Jij).

    Expected line format (whitespace-separated):
      i j rx ry rz Jij [dist]
    where i,j are integers (type indices) and Jij is a float (mRy).
    """
    pairs = []
    with open(jij_path, "r", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln:
                continue
            parts = ln.split()
            try:
                if len(parts) >= 6:
                    i = int(parts[0])
                    j = int(parts[1])
                    rx = float(parts[2])
                    ry = float(parts[3])
                    rz = float(parts[4])
                    Jij = float(parts[5])
                    pairs.append((i, j, rx, ry, rz, Jij))
            except ValueError:
                # skip header or malformed lines
                continue
    return pairs


def parse_dijso(dij_path: str):
    """Parse a `dijso.out`-style file into a list of (i,j,Dx,Dy,Dz).

    Expected line format (whitespace-separated):
      i j rx ry rz Dx Dy Dz [dist]
    where Dx,Dy,Dz are the D-vector components.
    """
    dm = []
    with open(dij_path, "r", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln:
                continue
            parts = ln.split()
            try:
                if len(parts) >= 8:
                    i = int(parts[0])
                    j = int(parts[1])
                    # skip rx,ry,rz at parts[2:5]
                    Dx = float(parts[5])
                    Dy = float(parts[6])
                    Dz = float(parts[7])
                    dm.append((i, j, Dx, Dy, Dz))
            except ValueError:
                continue
    return dm


def _write_dmfile(payload: Dict, path: str) -> str:
    """Write DM (Dzyaloshinskii–Moriya) vectors: i j Dx Dy Dz"""
    dm_list = payload.get("dm", [])
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# i j Dx Dy Dz\n")
        for entry in dm_list:
            if len(entry) >= 5:
                i, j, Dx, Dy, Dz = entry[:5]
                fh.write(f"{int(i)} {int(j)} {Dx: .8f} {Dy: .8f} {Dz: .8f}\n")
    return path


def export_from_files(lattice_nml: Optional[str] = None, input_nml: Optional[str] = None,
                      clust: Optional[str] = None, out_dir: str = "uppasd_out",
                      simid: str = "rs_export", ncell: Tuple[int, int, int] = (1, 1, 1),
                      BC: Tuple[str, str, str] = ("P", "P", "P"), Sym: int = 1,
                      jij_file: Optional[str] = None, dij_file: Optional[str] = None,
                      default_type: int = 1, default_moment: float = 1.0) -> List[str]:
    """
    High-level helper: parse provided RS-LMTO inputs and write UppASD-style
    output files into `out_dir`.

    Use either (lattice_nml + input_nml) or (clust + input_nml).
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Prefer clust+input, otherwise lattice+input. If input_nml is not
    # provided, try to parse lattice.nml alone (coords only) and apply
    # default types/moments.
    if clust and input_nml:
        struct = parse_rs_clust(clust, input_nml)
    elif lattice_nml and input_nml:
        struct = parse_rs_lmto_nml(lattice_nml, input_nml)
    elif lattice_nml and not input_nml:
        # parse lattice only for fractional coordinates
        import f90nml

        lat = f90nml.read(lattice_nml).get("lattice", None)
        # f90nml returns Cogroup-like objects; access safely
        def _safe_get(obj, key):
            if obj is None:
                return None
            getter = getattr(obj, "get", None)
            if callable(getter):
                return getter(key)
            try:
                return obj[key]
            except (KeyError, TypeError, IndexError):
                return None

        crd1 = _safe_get(lat, "crd(:, 1)") or []
        crd2 = _safe_get(lat, "crd(:, 2)") or []
        crd3 = _safe_get(lat, "crd(:, 3)") or []
        atoms = []
        n = len(crd1)
        for idx in range(n):
            atoms.append((str(default_type), [crd1[idx], crd2[idx], crd3[idx]]))
        a1 = _safe_get(lat, "a(:, 1)")
        a2 = _safe_get(lat, "a(:, 2)")
        a3 = _safe_get(lat, "a(:, 3)")
        struct = {"unit_cell": [a1, a2, a3], "atoms": atoms}
    else:
        raise ValueError("Provide either (clust + input_nml), (lattice_nml + input_nml), or (lattice_nml alone)")

    payload = build_uppasd_payload(struct, moments=None)

    # if no moments were found, set default_moment per atom
    if not any(payload.get("moments", [])):
        payload["moments"] = [float(default_moment)] * len(payload.get("atoms", []))

    # parse Jij/Dij files if provided and attach to payload
    if jij_file:
        jlist = parse_jijso(jij_file)
        # jlist entries are (i,j,rx,ry,rz,Jij) where i/j are type ids
        payload["jij"] = jlist
    if dij_file:
        dm_list = parse_dijso(dij_file)
        payload["dm"] = dm_list

    # write files in out_dir
    pos_path = os.path.join(out_dir, "posfile")
    mom_path = os.path.join(out_dir, "momfile")
    j_path = os.path.join(out_dir, "jfile")

    written = []
    written.append(_write_posfile(payload, pos_path))
    written.append(_write_momfile(payload, mom_path))
    # write jfile even if empty
    written.append(_write_jfile(payload, j_path))
    # write DM file if present
    dm_path = os.path.join(out_dir, "dmfile")
    if payload.get("dm"):
        written.append(_write_dmfile(payload, dm_path))

    # attempt to extract cell from payload lattice
    cell = payload.get("lattice")
    # ensure cell is a 3x3-like list; fall back to identity if missing
    try:
        if not cell or any(row is None for row in cell) or len(cell) < 3:
            cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    except Exception:
        cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    inpsd_path = write_inpsd_dat(
        out_dir,
        simid=simid,
        ncell=ncell,
        BC=BC,
        cell=cell,
        Sym=Sym,
        posfile="./posfile",
        momfile="./momfile",
        exchange="./jfile",
    )
    written.append(inpsd_path)
    return written


# Note: `export_from_files` (the UppASD directory-writing variant) is defined
# above. The older single-prefix writer has been removed to avoid confusion.


def _print_next_steps_help():
    print("RS-LMTO -> UppASD exporter (starter)")
    print("")
    print("What I need from you to finish the exporter:")
    print(" - Example RS-LMTO files: lattice.nml, input.nml, clust (if used), and any *_out.nml or report.out that contains magnetic data.")
    print(" - An example UppASD structure file or the exact format you want: ordering, header, units, whether spins/moments are included, and the Jij format (if desired).")
    print(" - Mapping rules for species->default moment if necessary, and whether indices should be 0- or 1-based in Jij files.")
    print("")
    print("Example usage (after you provide files):")
    print("  from source.export.rs_lmto_to_uppasd import export_from_files")
    print("  export_from_files('lattice.nml', 'input.nml', out_prefix='my_uppasd')")


if __name__ == '__main__':
    _print_next_steps_help()
