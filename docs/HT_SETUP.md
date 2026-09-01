# HT-safe setup

The high-throughput setup API validates the source namelists before invoking
the legacy Fortran kernel, runs that kernel in an isolated worker, validates
the generated files, and publishes them transactionally. The Python API is
the machine-facing interface; the legacy `-S` flag remains available for
interactive use.

```python
from pathlib import Path

from rsutil.core import SetupSpec
from rsutil.setup import SetupService, write_setup_report

source = Path("examples/bccFe")
destination = Path("runs/bccFe-setup")
spec = SetupSpec(seed=42, rec_method="block", rec_steps=20)

result = SetupService().run(
    source,
    destination,
    spec,
    collision_policy="error",
)
if not result.success:
    raise RuntimeError(result.to_dict())

write_setup_report(result, destination / "rsutil-setup-report.json")
```

For the familiar command-line workflow, run from a directory containing
`lattice.nml` and `input.nml`:

```bash
rsutil -S --validate --json-report rsutil-setup-report.json
```

The report is optional and does not change the ordinary human-readable
console output. Its top-level schemas are `rsutil-setup-report-v1` and
`rsutil-run-v1`; diagnostic `code` values are stable machine-facing values,
separate from human messages. Deterministic structure/setup hashes exclude
  timestamps, while the manifest records actual operation start/finish times
  and the worker/build environment separately. Stochastic provenance lists
  active operations and their seeds; a missing seed is only a warning when a
  stochastic distortion or SQS operation is actually requested.

## Lattice count validation

Setup validation follows the RS-LMTO `fable_v4` meanings of the lattice
counts; they are not interchangeable:

- `ntot` sizes the primitive-cell `crd`, `izp`, `no`, `iu`, and `ib` tables.
- `nrec` sizes `irec` and identifies recursion representatives.
- `nbas` is the charge-distribution/shell count and may differ from `ntot`.
- `nbulk_bulk` is a mode input used to derive bulk/reference counts; it is
  not treated as a universal array size.

For exchange/reference structures, `crd(:, site)` is the primitive-cell site
table, `izp(site)` is its one-based RS electronic type index, and
`atoms.label(type)` is the symbolic type label. RS-utils preserves this
source mapping; it does not infer or reclassify RS types from chemistry.

For surface inputs (`calctype='S'`), the validator checks the source-derived
relationships `nrec = nlay` and `ntype = nbulk_bulk + nlay` when those fields
are available. For layered/interface inputs (`calctype='L'`), it checks
`nrec = nlay` and `ntype = 2*nbulk_bulk + nlay`. Atom labels are checked
against `ntype`, not `nrec`, because impurity and surface workflows can have
frozen/reference labels in addition to recursion entries. Cluster-dependent
index bounds remain delegated to the protected RS-LMTO kernel because the
cluster does not exist at preflight time.
