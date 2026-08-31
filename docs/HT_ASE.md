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

The setup helper serializes the same `lattice.nml`/`input.nml` source contract
used by `SetupService`, then delegates to that service.  It does not invent
RS-LMTO element potential files; provide those as explicit runner resources
when using a real legacy setup kernel:

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
    collision_policy="error",
)
if not result.success:
    raise RuntimeError(result.to_dict())
```

The adapter accepts fractional `RSStructure` coordinates when converting to
ASE, and supports `angstrom`, `bohr`, `nm`, and `pm` length units.  Unknown
units fail explicitly instead of being silently interpreted as Å.
