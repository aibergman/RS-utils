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
- [x] **HT-05** — Transactional and idempotent setup
- [x] **HT-06** — Setup validation, provenance, and machine reports
- [x] **HT-07** — Canonical periodic exchange pair identity
- [x] **HT-08** — Extract `-J` behind the exchange service
- [x] **HT-09** — Correct exchange symmetry reduction
- [x] **HT-10** — Transactional exchange job writer and manifests
- [x] **HT-11** — Result parser foundation and SCF model
- [x] **HT-12** — Exchange, DMI, and tensor result parsing
- [x] **HT-13** — Transport and damping result contracts
- [x] **HT-14** — Versioned TB model and ML export
- [x] **HT-15** — ASE structure interoperability
- [x] **HT-16** — Deterministic campaign specification and generation
- [x] **HT-17** — Slurm array rendering, status, and collection
- [x] **HT-18** — HT readiness closeout and pilot evidence

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

### HT-05 decisions

- Setup-service publication now stages validated worker candidates and
  destination preimages on the destination filesystem, then atomically
  replaces the complete file set. A commit failure restores every file that
  was already replaced and removes newly created files/directories where
  possible.
- The library API exposes `collision_policy="error"` (default) and
  `collision_policy="replace"`. The strict default refuses existing output
  files; the legacy `-S` CLI explicitly selects `replace` for its familiar
  interactive behavior. In-place replacement of the two source namelists is
  authorized by `in_place=True`, while unrelated generated artifacts remain
  collision-protected unless replacement is explicitly selected.
- Backups are created only as part of the post-validation commit transaction,
  so staging, worker, validation, and commit failures leave source inputs
  byte-identical. Seed metadata records whether the setup is reproducible;
  any unspecified setup/SQS seed is marked non-reproducible.

### HT-06 decisions

- Added versioned `ValidationResult`, `RunManifest`, and `SetupReport` models
  with JSON-compatible serialization. Setup validation is exposed through the
  Python API and is applied both before the worker and to generated candidates.
- Semantic checks cover namelist structure, lattice geometry, counts/labels,
  recursion fields, supercell expansion, 2D output flags, generated output
  readability, and radius/sphere records. Diagnostics have stable codes and
  severity values; human console text is not used as a downstream contract.
- Added optional legacy CLI `--json-report PATH` output and documented the
  HT-safe Python/CLI setup usage in `docs/HT_SETUP.md`. Manifests record raw
  input/output hashes, deterministic structure/setup identity hashes, software
  and platform metadata, timestamps, validation evidence, and seed provenance.
- The legacy kernel currently emits `lld=21` for the default setup even when
  the normalized Python spec carries `rec_steps=20`. This remains a warning
  rather than a failure to preserve valid historical `-S` behavior; invalid
  recursion values and missing fields still fail semantically.

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

### HT-07 decisions

- Added immutable, zero-based `PairKey(i, j, R)` and explicit one-based
  RS/Fortran conversion helpers. `R` is always an integer translation in the
  reference lattice basis; Cartesian geometry is derived separately.
- Added `ExchangePair` and a dependency-light mapper that resolves cluster
  centers/neighbors against `RSStructure` sites by species and lattice
  geometry. Non-unique or inconsistent matches fail explicitly, and mapped
  pairs are sorted by canonical key.
- Preserved the legacy `R_ij` mapping output, but label it explicitly as a
  Cartesian displacement. No symmetry reduction or `-J` generation behavior
  was changed; those remain HT-08/HT-09 scope.

### HT-08 decisions

- Added a path-explicit exchange service with focused reference-input,
  structured Fermi parsing, neighbor discovery, canonical job planning, and
  compatibility serialization modules. `ExchangePlan` exposes the canonical
  `PairKey` inventory while `ExchangeJobSpec` retains the legacy one-based
  cluster indices needed in generated `Atom_*` namelists.
- Legacy `exchange_map()` and helper names remain as thin wrappers. The `-J`
  CLI now routes through `ExchangeService`, honors `--data` as its source
  directory, and reports exchange/parser/output failures as nonzero CLI
  operations. Explicit mapping output is attempted before job publication so
  a requested mapping failure does not create a successful-looking job set.
- The pre-HT-09 `--unique`/`--soc` reduction algorithm remains available
  behind an isolated reducer hook. Its symmetry-operation filtering is now
  computed from public pymatgen operation data rather than mutating the
  analyzer's private state; scientific reduction correctness remains HT-09
  scope.

### HT-09 decisions

- Replaced shell-width/vector heuristics with exact periodic pair-key orbits.
  Public pymatgen symmetry operations are compiled into site permutations and
  integer cell-image shifts, applying `(i,j,R) -> (p(i),p(j),n(j)+W R-n(i))`.
- Legacy `--soc` now means a collinear `+z` Cartesian magnetization reference.
  Operations are filtered using axial-vector transformation and species/site
  preservation; generic noncollinear magnetic symmetry is unsupported and is
  not inferred from the Boolean flag.
- Added `--full-pairs`/`full_pairs=True` as a validation opt-out. The legacy
  `--shell-width` parameter remains accepted for compatibility but has no role
  in scientific equivalence. See `docs/HT_EXCHANGE_SYMMETRY.md`.

### HT-10 decisions

- Extended `ExchangeJobSpec` with the canonical pair inventory, cutoff and
  reduction/SOC policy, explicit Fermi-energy source, required converged atom
  files, and source-run identity/hash. Required `<label>_out.nml` files are
  derived from the selected input namelist and copied by explicit name from
  the selected source directory; unrestricted CWD globs are no longer used.
- Exchange jobs are fully staged and semantically validated before publication.
  The versioned `rsutil_exchange_manifest.json` records the canonical pair
  list, one-based legacy `ijpair` semantics, deterministic specification hash,
  source-run hash, and generated-file hashes.
- The library exposes `error`, `replace_generated`, and `reuse_if_identical`
  collision policies. Only valid rsutil-owned directories may be replaced or
  reused by default; unmarked legacy/user directories require explicit
  `force_legacy_collision=True`. Legacy `-J` retains the `Atom_*` layout and
  uses the owned-directory replacement policy. See `docs/HT_EXCHANGE_JOBS.md`.

### HT-11 decisions

- Added the dependency-light `rsutil.results` package with versioned
  `ElectronicResult`/`SCFResult` models, explicit `Quantity` units,
  `ParserProvenance` source hashes, separate parse diagnostics and validation,
  and deterministic JSON-compatible serialization.
- `parse_run(path)` accepts either an explicit run directory or `report.out`.
  It parses the stable report records for Fermi energy, energies, moments,
  orientations, occupations, charge transfer, and RMS difference while
  retaining unknown units where the upstream text contract does not establish
  a unit. Atom labels are converted from the report's one-based form to the
  core's zero-based site index.
- Convergence is populated only from explicit `Converged!` or `Not converged!`
  markers in `out.scf`/`out.cold`; a report alone yields `unknown`. The legacy
  exchange Fermi reader is now a compatibility wrapper around this parser.
- Material evidence limitation: this checkout contains the checked-in
  B2Fe/FePt reports and a local `fable_v3` source checkout, while no fable_v4
  source is present in the repository. The inspected report-writing contract
  does not emit convergence or units for total energies/moments, so those
  fields remain explicitly unknown rather than being inferred.

### HT-12 decisions

- Added versioned `InteractionRecord`/`ExchangeSet` models keyed by canonical `PairKey(i,j,R)`, with independent optional isotropic exchange, DMI-vector, and full-tensor quantities. JSON export preserves units, source rows/hashes, channel metadata, and optional derived pair geometry.
- Added strict parsers for the real checked-in `jijso.out`/`jijfo.out` scalar and `dijso.out`/`dijfo.out` DMI table schemas. Rows are mapped by one-based output indices plus Cartesian displacement against the exchange manifest/reference structure; raw row order and `Atom_*` names are not used as identity. Ambiguous or mismatched mapping fails explicitly.
- Added explicit tensor-shaped row support only when a filename/header identifies a tensor, plus separate reverse-pair physics validation (J symmetric, D antisymmetric, tensor transpose) that reports warnings by default and does not impose crystal constraints.
- Material evidence limitation: this checkout contains no `fable_v4/source/exchange.f90` or real fable_v4 tensor fixture. The scalar/DMI documentation is based on the checked-in MnSi outputs and `docs/notes/RSLMTO_format.md`; tensor parsing remains fail-closed pending upstream fable_v4 evidence.

### HT-13 decisions

- The checked-out RS-utils tree contains no `fable_v4` source, conductivity/Kubo/KPM implementation, Gilbert-damping implementation, transport/damping writer, or corresponding fixture. No transport or damping text format is therefore production-supported; the evidence and future upstream requirement are recorded in `docs/HT_TRANSPORT_DAMPING_FORMATS.md`.
- Added `rsutil.results.transport` with explicit `ConductivityResult` and `DampingResult` contracts. They require units, axes, component/index ordering, method, channel, approximation, and parser provenance, and validate the declared value shape.
- Added a deliberately small registry and fail-closed parser entry points. Candidate files return `UNSUPPORTED_AMBIGUOUS` with output/input/lattice/manifest hashes and no guessed values; missing files remain explicit missing-file errors. No RS-LMTO scientific source was changed.

### HT-14 decisions

- Verified the external `fable_v4/source/hamiltonian_paoflow_io.f90` writer and
  implemented its `rs-tb-meta-v1`/`rs-tb-hoppings-v1` pair. Metadata uses
  Cartesian lattice/site coordinates, eV hopping values, one-based text
  indices, and hopping records keyed by integer target-cell translation plus
  site/local/global orbital fields. The exact schema is documented in
  `docs/HT_TIGHT_BINDING.md`.
- Added a dependency-light `TightBindingModel` with deterministic `H(R)`
  identity/order, complex `H(k)` reconstruction, bands, provenance, and
  optional Fermi metadata. The upstream versioned writer does not emit an
  overlap matrix, so overlap remains absent rather than guessed.
- Added deterministic compressed NPZ plus JSON ML export preserving
  `(R, site_i, site_j, orbital_i, orbital_j)` and native complex values. No
  HDF5/Zarr dependency was added. Legacy `rs2txt.dat` reading remains
  unchanged; conversion is available only when periodic identity is
  unambiguous.
- `control%nsp` is not present in the current TB metadata writer. When an
  input namelist is supplied, RS-utils records the explicit RS-LMTO rule:
  `nsp=1` is magnetic collinear without SOC (GBT-compatible), while `nsp=2`
  or `4` enables SOC. Missing `nsp` is reported as not encoded.
- The upstream repository provides the writer but no committed generated
  versioned TB output fixture; the local fixture is an exact-schema,
  complex-valued regression fixture derived from that writer contract.

### HT-15 decisions

- Added an isolated, lazy-importing `rsutil.ase` adapter.  ASE Cartesian
  positions/cells are recorded in explicit Å units; cell rows, PBC flags,
  species, atom order, JSON-compatible `Atoms.info`, and per-atom arrays are
  preserved without placing ASE objects in core contracts or manifests.
- Added `ASESetupSource`/`run_setup_from_ase`, which serializes only the
  existing `lattice.nml`/`input.nml` source contract and delegates execution to
  `SetupService`.  Element potential files are not guessed and must be passed
  as explicit runner resources.  No ASE calculator was added.

### HT-16 decisions

- Added the versioned JSON `rsutil-campaign-v1` specification and typed
  `Campaign`, `CampaignEntry`, and `CampaignIndex` contracts.  Entry IDs are
  SHA-256 hashes of source file contents plus normalized setup, calculation,
  exchange, and entry metadata; entry order and human names do not determine
  scientific identity.
- Campaign preflight validates every source, setup specification, namelist,
  and required stochastic seed before creating the destination.  Generation
  composes `SetupService` once per entry and writes relative-path
  `run_manifest.json` files and a relocatable `rsutil-campaign-index-v1`
  index.  Existing identical runs are reused by manifest hash; incompatible
  or unmarked directories are explicit collisions.
- The primary CLI surface is `rsutil campaign validate SPEC` and `rsutil
  campaign create SPEC [--destination PATH]`.  YAML is accepted only as a
  convenience when PyYAML is already installed; no scheduler or job
  execution is introduced.  Invalid preflight entries leave the campaign
  destination untouched; runtime failures retain per-entry setup atomicity.

### HT-17 decisions

- Added `rsutil.campaign.render_slurm_array`, which consumes the relocatable
  HT-16 index and writes a stable zero-based TSV array map plus a readable,
  shell-quoted `sbatch` script. Executable, launcher, environment preamble,
  partition/account/time, and resource directives are all render parameters;
  no site-specific scheduler or module behavior is selected.
- Added filesystem-only `status_campaign` classification. Explicit failure
  and lock markers are recognized, parser/validation evidence classifies
  completed runs, and missing/invalid evidence remains conservative. Slurm is
  not queried, so the report explicitly states `scheduler_state=not_queried`.
- Added failure-tolerant `collect_campaign` aggregation with versioned JSON,
  optional CSV/JSONL scalar summaries, structured result sidecars, raw-file
  hashes, manifest identity/content hashes, and exchange-subjob canonical pair
  references. No retry, restart, dependency, daemon, or workflow-engine
  behavior was added; those remain a future AiiDA boundary.

### HT-18 decisions

- Closed the foundation with the measured readiness report in
  `docs/HT_READINESS_REPORT.md` and the reproducible four-entry pilot
  specification in `examples/campaigns/ht18_pilot.json`.
- The pilot evidence combines real setup generation, checked-in electronic,
  exchange/DMI, versioned TB, and fail-closed transport fixtures. A stable
  transport/damping output format is not present in this checkout, so that
  mode remains explicitly unsupported rather than being promoted by the
  closeout.
- The report records the direct-vs-worker golden comparison, legacy `-S`/`-J`
  semantic corpus, atomicity/failure-injection tests, deterministic campaign
  regeneration, and fake-executable Slurm/status/collection dry-run. AiiDA
  and UppASD workflow integration remain outside scope.
