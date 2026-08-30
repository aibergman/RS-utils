# RS-utils High-Throughput Foundation Blueprint

This is the repository master ledger for the RS-utils high-throughput (HT)
campaign. It is governed by
`docs/dev/RSUTIL_HT_LUNA_BLUEPRINT_PACK/00_HEAD_RSUTIL_HT_FOUNDATION.md`.

## Mission and scope

Refactor RS-utils from a collection of interactive utilities into a reliable
RS-LMTO high-throughput infrastructure layer while preserving the useful
legacy workflows. The campaign targets RS-LMTO electronic, magnetic,
exchange, transport, and real-space tight-binding data, with typed contracts,
isolated legacy-kernel execution, deterministic job generation, provenance,
and lightweight campaign tooling.

The central design rule is **aggressive internal hardening with deliberately
conservative external compatibility**.

The external RS-LMTO reference for this campaign is the `fable_v4` branch.
That repository is not modified by this task.

## Compatibility policy

- **BC-01:** scientifically valid `rsutil -S ...` and `rsutil -J ...`
  invocations remain supported with familiar defaults and ordinary file and
  directory conventions.
- **BC-02:** `-S`/`-J` remain available; a new subcommand interface is not a
  prerequisite for HT work.
- Compatibility is semantic and workflow-oriented, not byte-for-byte output
  identity.
- Existing band/DOS, visualization, LMTO47, and UppASD-export utilities are
  outside the refactoring scope except where shared infrastructure requires
  smoke coverage.

The detailed contract is in [`HT_COMPATIBILITY_CONTRACT.md`](HT_COMPATIBILITY_CONTRACT.md).
The checked-out-code inventory is in [`HT_BASELINE_INVENTORY.md`](HT_BASELINE_INVENTORY.md).

## Protected Fortran-kernel policy

The legacy Fortran sphere, empty-sphere, lattice/input I/O, and potential
parameter routines are a protected numerical compatibility kernel. Their
algorithms, tolerances, and formulas are not ported or modernized as part of
this campaign. HT safety is achieved by containment: future setup services
invoke the existing `ctypes` entry point in an isolated worker workspace so
`STOP`/`ERROR STOP` paths cannot terminate the parent process, and fixed
filenames remain local to that workspace.

## Task tracking

- [x] **HT-00** — Freeze baseline and compatibility contract
- [x] **HT-01** — Legacy semantic regression harness
- [x] **HT-02** — Core contracts and safe path primitives
- [x] **HT-03** — Isolate the legacy Fortran setup kernel
- [x] **HT-04** — Extract `-S` behind the setup service
- [ ] **HT-05** — Transactional and idempotent setup
- [ ] **HT-06** — Setup validation, provenance, and machine reports
- [ ] **HT-07** — Canonical periodic exchange pair identity
- [ ] **HT-08** — Extract `-J` behind the exchange service
- [ ] **HT-09** — Correct exchange symmetry reduction
- [ ] **HT-10** — Transactional exchange job writer and manifests
- [ ] **HT-11** — Result parser foundation and SCF model
- [ ] **HT-12** — Exchange, DMI, and tensor result parsing
- [ ] **HT-13** — Transport and damping result contracts
- [ ] **HT-14** — Versioned TB model and ML export
- [ ] **HT-15** — ASE structure interoperability
- [ ] **HT-16** — Deterministic campaign specification and generation
- [ ] **HT-17** — Slurm array rendering, status, and collection
- [ ] **HT-18** — HT readiness closeout and pilot evidence

## Dependency summary

```text
HT-00 -> HT-01 -> HT-02 -> HT-03 -> HT-04 -> HT-05 -> HT-06
                                                        |
                 +------------------+-------------------+------------------+
                 |                  |                                      |
                 v                  v                                      v
               HT-07              HT-11                                  HT-15
                 |                  |                                      |
                 v             +----+----+                                 |
               HT-08           v         v                                 |
                 |            HT-13     HT-14                              |
                 v              \         /                                 |
               HT-09 -> HT-10 ---+-------+---------------------------------+
                                      |
                                      v
                              HT-16 -> HT-17 -> HT-18
```

HT-12 consumes the pair/job-manifest work from HT-10 and parser contracts
from HT-11. HT-13 follows HT-11. HT-14 follows the parser/data contracts and
the `fable_v4` TB reference. HT-16 composes the earlier setup, exchange,
parser, and provenance APIs; it is not a general workflow engine.

## Scope boundaries and non-goals

- No AiiDA plugin, new UppASD workflow redesign, or general workflow engine.
- No Python port or broad rewrite of the Fortran sphere/potential-parameter
  kernel.
- No ML model training or database server/backend project.
- No gratuitous CLI redesign and no unrelated style refactor.
- No RS-LMTO scientific-kernel changes merely to make parsing easier; missing
  stable upstream output should be documented as an external requirement.

## Decisions and deviations

### HT-00 baseline decisions

- The source checkout reports version `0.10.0`, and the baseline commit is
  recorded in `HT_BASELINE_INVENTORY.md`.
- `rsutil --help` resolves to an external installed launcher in this
  environment and fails before reaching this checkout. The equivalent source
  entry point, `PYTHONPATH=. python -m rsutil.cli --help`, was used and its
  relevant output is archived in the inventory.
- The source checkout successfully loads its built shared library through
  `ctypes`; its observed build path is recorded in the inventory.
- A direct source-checkout `-J` probe reaches a current parser defect: the
  code reads `args.rij_out` and `args.rij_format`, but the parser does not
  define those options. This is recorded for a later focused compatibility
  fix and is intentionally not changed by HT-00.
- The checkout has no collected pytest tests. Existing CI builds the project;
  the pip workflow performs a package/example smoke invocation rather than a
  pytest suite. Historical claims in `docs/notes/TEST_RESULTS.md` are
  retained as notes, not treated as executable baseline evidence.

No scientific or functional code was changed for HT-00.

### HT-02 decisions

- Added a dependency-light `rsutil.core` package with library-native
  `RSStructure`/site contracts, `SetupSpec`, structured diagnostics and
  operation results. Site indices are canonical zero-based contiguous values;
  coordinate convention and length unit are explicit.
- Scientific setup parameters are kept separate from output and in-place
  policy. Deterministic JSON/specification hashes exclude temporal provenance
  fields, while file hashes cover content only.
- Added owned temporary/staging workspaces, explicit file staging, and atomic
  replacement helpers. They never change process CWD; same-filesystem staging
  is created beside the eventual destination.
- Made the existing top-level `bands` convenience import lazy so the
  dependency-light core can be imported without loading optional plotting
  dependencies. The public `rsutil.bands` attribute remains available.

### HT-03 decisions

- Kept the numerical Fortran sources unchanged. The direct ctypes conversion
  is now explicitly named `_process_lattice_unsafe`; `rsutil.worker` invokes
  it only in a child process after changing into an owned workspace.
- Added `FortranSetupRunner`, which stages only the two files used by the
  inspected setup call graph (`lattice.nml` and `input.nml`), accepts explicit
  future resources, captures worker output/exit status, and inventories
  generated namelists and sphere files. It returns `SetupResult` diagnostics
  and does not commit outputs; setup-service commit policy remains HT-04/05.
- The runner defaults to no timeout so valid scientific setup work is not
  rejected arbitrarily; callers can provide a positive timeout. A caller
  supplied workspace is retained for inspection/commit, while auto-created
  workspaces are cleaned up.

### HT-04 decisions

- Added `SetupService`/`run_setup` with explicit source and destination paths;
  it owns setup-spec normalization, validation, supercell/SQS transforms,
  worker orchestration, output inventory, and candidate-file publication.
- Legacy `-S` now translates its existing flags into `SetupSpec`, calls the
  service, and retains the expected `lattice_new.nml`, `input_new.nml`, sphere
  outputs, console presentation, seeded SQS behavior, and `--in-place` flow.
- Setup transforms and the safe worker never change the reusable caller's
  process CWD. The legacy `process_lattice()` wrapper remains an explicit
  compatibility entry point for direct in-process callers.
- Material deviation: setup candidate publication is implemented in HT-04;
  the stronger rollback/idempotence guarantees for partial commit failures
  remain the focused scope of HT-05.

### HT-01 regression-harness decisions

- The pytest suite stages only the selected files from the existing bccFe,
  FeN, and FePt examples into temporary directories; tests never modify the
  repository examples in place.
- The baseline `-J` parser referenced undocumented-in-code `--rij-out` and
  `--rij-format` attributes. These options were restored so the documented
  mapping workflow is executable, and explicitly requested mapping-write
  failures now produce a failed operation instead of a printed success path.
- `--unique`/`--soc` coverage freezes availability and non-crash behavior only.
  Scientific pair-equivalence expectations remain deferred to HT-09.
