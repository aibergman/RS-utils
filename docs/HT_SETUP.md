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
timestamps, while the manifest records timestamps and the worker/build
environment separately.
