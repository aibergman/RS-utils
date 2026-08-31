# HT-15 — ASE interoperability

ASE is an optional boundary adapter.  `RSStructure` remains the native
RS-utils contract, and ASE objects are not stored in manifests or result
models.  Install the optional dependency with:

```bash
pip install 'rsutil[ase]'
```

Conversion preserves cell-row orientation, Cartesian positions in Å, PBC
flags, species, atom order, and JSON-compatible ASE metadata/arrays:

```python
from ase.io import read

from rsutil.ase import from_ase, to_ase

atoms = read("structure.cif")
structure = from_ase(atoms)
round_tripped = to_ase(structure)
```

The `Atoms <-> RSStructure` conversion above is the supported HT ASE
boundary. Full ASE-to-setup-source generation is retained only as an
explicitly opt-in experimental convenience; it emits placeholder LMTO policy
(`rc`, `ct`, recursion dimensions, and related values) and does not establish
scientifically reviewed setup defaults. It also does not invent RS-LMTO
element potential files. If this experimental path is deliberately reviewed,
pass `experimental=True` and provide the required resources explicitly:

```python
from ase.io import read

from rsutil.ase import run_setup_from_ase
from rsutil.core import SetupSpec

atoms = read("structure.cif")
result = run_setup_from_ase(
    atoms,
    "runs/structure-setup",
    SetupSpec(seed=42),
    resources=["potentials/Fe.nml", "potentials/Pt.nml"],
    experimental=True,
    collision_policy="error",
)
if not result.success:
    raise RuntimeError(result.to_dict())
```

The adapter accepts fractional `RSStructure` coordinates when converting to
ASE, and supports `angstrom`, `bohr`, `nm`, and `pm` length units.  Unknown
units fail explicitly instead of being silently interpreted as Å.

The experimental source helper derives the initial Wigner-Seitz radius as
`(3 V_atom / (4 pi))^(1/3)`, but this value and the other generated policy
fields remain subject to review by the legacy setup workflow.
