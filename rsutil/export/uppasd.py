"""
RS-LMTO -> UppASD exporter (moved to rsutil.export.uppasd)

This is the same starter exporter previously located at
`source/export/rs_lmto_to_uppasd.py` but moved into the package so it can be
imported as `rsutil.export.uppasd` and exposed from the package CLI.
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
    import f90nml

    if InputReader is not None and hasattr(InputReader, "read_rs_lmto_nml"):
        struct = InputReader.read_rs_lmto_nml(lattice_nml, input_nml)
        return struct

    lat = f90nml.read(lattice_nml).get("lattice", {})
    inp = f90nml.read(input_nml)

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
    while True:
        key = f"label({i})"
        if key in atoms_block:
            raw = atoms_block[key]
            species.append("".join(ch for ch in str(raw).strip('\'\"') if ch.isalpha()))
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
    if InputReader is not None and hasattr(InputReader, "read_rs_clust"):
        return InputReader.read_rs_clust(clust_file, input_file)

    atoms = []
    with open(clust_file, "r", encoding="utf-8") as f:
        lines = [ln.rstrip() for ln in f if ln.strip()]
    if not lines:
        raise ValueError("Empty clust file")
    import f90nml

    inp = f90nml.read(input_file)
    atoms_block = inp.get("atoms", {})
    species_lookup = []
    i = 1
    while True:
        key = f"label({i})"
        if key in atoms_block:
            raw = atoms_block[key]
            species_lookup.append("".join(ch for ch in str(raw).strip('\'\"') if ch.isalpha()))
            i += 1
        else:
            break

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
    payload = {"lattice": None, "atoms": [], "moments": [], "jij": []}
    if Structure is not None and isinstance(structure, Structure):
        payload["lattice"] = getattr(structure, "unit_cell", None)
        atoms_obj = getattr(structure, "atoms_uc", None)
        if atoms_obj is None:
            atoms_obj = getattr(structure, "atoms", None)
        if atoms_obj:
            for atom in atoms_obj:
                spec = getattr(atom, "species", getattr(atom, "atom_type", "X"))
                pos = getattr(atom, "position", None)
                payload["atoms"].append((spec, list(pos)))
    else:
        payload["lattice"] = structure.get("unit_cell") if isinstance(structure, dict) else None
        payload["atoms"] = structure.get("atoms", []) if isinstance(structure, dict) else []

    if moments is not None:
        payload["moments"] = list(moments)
    else:
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


def _write_posfile(payload: Dict, path: str) -> str:
    atoms = payload.get("atoms", [])
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
    atoms = payload.get("atoms", [])
    moments = payload.get("moments", [0.0] * len(atoms))
    with open(path, "w", encoding="utf-8") as fh:
        for idx, m in enumerate(moments, start=1):
            chem = 1
            ux, uy, uz = 1.0, 0.0, 0.0
            fh.write(f"{idx} {chem} {m: .6f} {ux: .6f} {uy: .6f} {uz: .6f}\n")
    return path


def _write_jfile(payload: Dict, path: str) -> str:
    jlist = payload.get("jij", [])
    with open(path, "w", encoding="utf-8") as fh:
        for entry in jlist:
            if len(entry) == 3:
                i, j, Jij = entry
                rx, ry, rz = 0.0, 0.0, 0.0
            elif len(entry) >= 6:
                i, j, rx, ry, rz, Jij = entry[:6]
            else:
                continue
            fh.write(f"{int(i)} {int(j)}  {rx: .6f}  {ry: .6f}  {rz: .6f}  {Jij: .14f}\n")
    return path


def write_inpsd_dat(out_dir: str, simid: str = "rs_export", ncell: Tuple[int, int, int] = (1, 1, 1),
                    BC: Tuple[str, str, str] = ("P", "P", "P"), cell: Optional[object] = None,
                    Sym: int = 1, posfile: str = "./posfile", momfile: str = "./momfile",
                    exchange: str = "./jfile", dm: Optional[str] = None) -> str:
    path = os.path.join(out_dir, "inpsd.dat")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"simid {simid}\n")
        fh.write(f"ncell     {ncell[0]:d}        {ncell[1]:d}       {ncell[2]:d}\n")
        fh.write(f"BC        {BC[0]}         {BC[1]}         {BC[2]}\n")
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
                continue
    return pairs


def parse_dijso(dij_path: str):
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
                    Dx = float(parts[5])
                    Dy = float(parts[6])
                    Dz = float(parts[7])
                    dm.append((i, j, Dx, Dy, Dz))
            except ValueError:
                continue
    return dm


def _write_dmfile(payload: Dict, path: str) -> str:
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
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if clust and input_nml:
        struct = parse_rs_clust(clust, input_nml)
    elif lattice_nml and input_nml:
        struct = parse_rs_lmto_nml(lattice_nml, input_nml)
    elif lattice_nml and not input_nml:
        import f90nml

        lat = f90nml.read(lattice_nml).get("lattice", None)

        # f90nml typically stores 'a' as a list-of-lists and 'crd' as a list of
        # coordinate triples. Handle both that format and the older 'crd(:, 1)'
        # style.
        a = None
        crd = None
        if lat is not None:
            # prefer the simple keys
            a = getattr(lat, 'get', lambda k, d=None: None)('a') or lat.get('a') if hasattr(lat, 'get') else None
            crd = getattr(lat, 'get', lambda k, d=None: None)('crd') or lat.get('crd') if hasattr(lat, 'get') else None
            # fallback to explicit names used in some files
            if crd is None:
                # try 'crd(:, 1)' format
                try:
                    c1 = lat.get("crd(:, 1)")
                    c2 = lat.get("crd(:, 2)")
                    c3 = lat.get("crd(:, 3)")
                    if c1 and c2 and c3:
                        crd = [[c1[i], c2[i], c3[i]] for i in range(len(c1))]
                except Exception:
                    crd = None
            if a is None:
                try:
                    a1 = lat.get("a(:, 1)")
                    a2 = lat.get("a(:, 2)")
                    a3 = lat.get("a(:, 3)")
                    if a1 and a2 and a3:
                        a = [a1, a2, a3]
                except Exception:
                    a = None

        atoms = []
        if crd:
            # If crd entries are fractional coordinates and a is provided
            # convert to Cartesian: r_cart = frac * A (A rows are lattice vectors)
            if a and isinstance(a, list) and len(a) >= 3:
                A = a
                for pt in crd:
                    # pt is [u,v,w]
                    x = pt[0] * A[0][0] + pt[1] * A[1][0] + pt[2] * A[2][0]
                    y = pt[0] * A[0][1] + pt[1] * A[1][1] + pt[2] * A[2][1]
                    z = pt[0] * A[0][2] + pt[1] * A[1][2] + pt[2] * A[2][2]
                    atoms.append((str(default_type), [x, y, z]))
            else:
                for pt in crd:
                    atoms.append((str(default_type), [pt[0], pt[1], pt[2]]))

        struct = {"unit_cell": a, "atoms": atoms}
    else:
        raise ValueError("Provide either (clust + input_nml), (lattice_nml + input_nml), or (lattice_nml alone)")

    payload = build_uppasd_payload(struct, moments=None)

    if not any(payload.get("moments", [])):
        payload["moments"] = [float(default_moment)] * len(payload.get("atoms", []))

    if jij_file:
        jlist = parse_jijso(jij_file)
        payload["jij"] = jlist
    if dij_file:
        dm_list = parse_dijso(dij_file)
        payload["dm"] = dm_list

    pos_path = os.path.join(out_dir, "posfile")
    mom_path = os.path.join(out_dir, "momfile")
    j_path = os.path.join(out_dir, "jfile")

    written = []
    written.append(_write_posfile(payload, pos_path))
    written.append(_write_momfile(payload, mom_path))
    written.append(_write_jfile(payload, j_path))
    dm_path = os.path.join(out_dir, "dmfile")
    if payload.get("dm"):
        written.append(_write_dmfile(payload, dm_path))

    cell = payload.get("lattice")
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


def _print_next_steps_help():
    print("RS-LMTO -> UppASD exporter (starter)")
    print("")
    print("What I need from you to finish the exporter:")
    print(" - Example RS-LMTO files: lattice.nml, input.nml, clust (if used), and any *_out.nml or report.out that contains magnetic data.")
    print(" - An example UppASD structure file or the exact format you want: ordering, header, units, whether spins/moments are included, and the Jij format (if desired).")
    print(" - Mapping rules for species->default moment if necessary, and whether indices should be 0- or 1-based in Jij files.")
    print("")
    print("Example usage (after you provide files):")
    print("  from rsutil.export.uppasd import export_from_files")
    print("  export_from_files('lattice.nml', 'input.nml', out_dir='my_uppasd')")


if __name__ == '__main__':
    _print_next_steps_help()
