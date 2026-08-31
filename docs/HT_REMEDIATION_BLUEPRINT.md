# RS-utils HT remediation blueprint

## Status

**REMEDIATION OPEN — the current branch is not HT-ready for scientific
production.**

This ledger reopens the historical HT-18 closeout after the adversarial audit.
The existing HT architecture and HT-18 evidence remain useful baseline
material, but they do not close the release gates below. Later remediation
tasks must preserve the findings and acceptance criteria recorded here.

The external scientific reference remains the public RS-LMTO-ASA `fable_v4`
branch. RS-LMTO `fable_v4` is not modified by this campaign.

## Governing invariants

### FORTRAN-R1 — Protected numerical kernel

Do not modify numerical algorithms, constants, tolerances, or formulas in the
legacy Fortran sphere/potential-parameter kernel during remediation.

The current subprocess isolation architecture should remain.

### BC-R1 — Legacy command compatibility

Scientifically valid uses of:

```bash
rsutil -S ...
rsutil -J ...
```

must continue to work with familiar CLI behavior and ordinary file/directory
conventions, except where the previous behavior was scientifically incorrect
or unsafe.

### DATA-R1 — Explicit source units

Parsed values must carry the unit established by the source file contract.
Never silently rescale a source value merely to make internal conventions
convenient. If canonical converted values are offered, preserve source
value/unit provenance.

### STATE-R1 — Do not collapse independent states

Execution success, parse validity, and scientific convergence are distinct.

### SYM-R1 — Reduction must be conservative

If symmetry equivalence cannot be proven for the supplied crystallographic
**and magnetic** reference state, retain the calculation rather than
eliminate it.

### PROV-R1 — Reuse means current bytes match scientific identity

A manifest hash alone is insufficient if the materialized scientific inputs
have been changed after generation.

### TEST-R1 — Scientific negative controls are mandatory

Every repaired bug must gain a regression test that would fail with the
pre-remediation implementation.

## Task tracking

- [x] **HT-R00** — Reopen readiness and freeze adversarial acceptance criteria
- [x] **HT-R01** — Correct exchange/DMI/tensor units and unit evidence
- [x] **HT-R02** — Restore complete RS-LMTO setup/lattice validation semantics
- [x] **HT-R03** — Implement scientific campaign state and execution evidence
- [x] **HT-R04** — Harden manifests, pair identity, and reuse integrity
- [x] **HT-R05** — Bound and harden exchange symmetry/magnetic-reference semantics
- [x] **HT-R06** — Complete transport/damping source audit against `fable_v4`
- [x] **HT-R07** — Provenance, ASE, path, indexing, and dead-code cleanup
- [x] **HT-R08** — Make packaging, test environment, and CI reproducible
- [ ] **HT-R09** — True adversarial HT closeout and pilot evidence *(superseded/reopened by post-R09 self-review; historical evidence preserved below)*
- [x] **HT-R10** — Interaction scientific contract closeout
- [x] **HT-R11** — Authoritative regeneration-based reuse
- [x] **HT-R12** — Magnetic-reference-aware symmetry contract
- [x] **HT-R13** — Installed-wheel clean-environment verification
- [x] **HT-R14** — Final falsification closeout

No task may be marked complete before its acceptance tests pass. The historical
R09 checkbox is reopened because the post-R09 self-review found additional
release blockers; its evidence below is retained for traceability only.

## Frozen adversarial blockers

| Severity | Affected module | Observed failure | Remediation task | Release gate |
| --- | --- | --- | --- | --- |
| P0-A | `rsutil/results/units.py`, `rsutil/results/exchange.py` | Scalar exchange/DMI values from RS-LMTO are labeled `Ry` although the writer emits them in `mRy`; the real fixture test preserves the factor-of-1000 error. | HT-R01 | Exchange source units are correct and supported by current `fable_v4` source/writer evidence. |
| P0-B | `rsutil/setup_validation.py` | `nbulk_bulk`, `ntot`, `nbas`, and `nrec` are treated as universally equal to `nrec`, rejecting valid hcp, surface, impurity, defect, and other modes. | HT-R02 | Setup validation follows source-defined RS-LMTO semantics and accepts every claimed supported mode. |
| P0-C | `rsutil/campaign/ht17.py`, result/status models | Parser/validation success can produce `COMPLETED_VALID` even when `report.out` explicitly says the SCF did not converge. | HT-R03 | Execution outcome, parse validity, and scientific convergence are represented independently; non-converged runs are never scientifically successful. |
| P0-C | Generated campaign runner and persisted run evidence | Runner execution exit status is not robustly persisted, so a failed execution can lack durable failure evidence. | HT-R03 | Failed execution always leaves persistent, machine-readable exit evidence. |
| P1 | `rsutil/exchange/jobs.py`, exchange manifests and reuse path | Pair identity has duplicated representations, and `reuse_if_identical` does not prove that materialized scientific inputs still match the recorded identity after edits. | HT-R04 | Pair identity is authenticated and reuse verifies current scientific input bytes. |
| P1 | `rsutil/exchange/symmetry.py` | Reduction is more general in presentation than in proven scope and lacks a complete site-resolved magnetic-state contract for AFM, ferrimagnetic, and noncollinear references. | HT-R05 | Unproven magnetic equivalence retains jobs; only the documented proven domain is reduced. |
| P1 | `rsutil/results/transport.py`, `docs/HT_TRANSPORT_DAMPING_FORMATS.md` | Transport/damping support was closed conservatively without the required current `fable_v4` source audit; no stable machine-readable writer contract is evidenced. | HT-R06 | Transport/damping domains are source-audited and explicitly classified as ready, restricted, experimental, or unsupported. |
| P1 | `rsutil/setup_service.py`, `rsutil/setup_validation.py`, campaign generation | Deterministic setups with no stochastic operation can be marked non-reproducible merely because no seed is supplied. | HT-R07 | Reproducibility provenance reflects actual stochastic operations and supplied seeds. |
| P1 | `rsutil/ase.py` setup helper | ASE-to-full-setup helper invents LMTO policy and computes `wav` as `V_atom^(1/3)`, not the Wigner-Seitz-radius formula. | HT-R07 | ASE structure interop remains useful, while setup policy is explicit, reviewed, and scientifically correct or quarantined. |
| P1 | `pyproject.toml`, `.github/workflows/`, test environment | Clean-checkout evidence is not fully reproducible: dependencies are not cleanly locked/grouped and the pytest CI path does not run on all HT-branch pushes. | HT-R08 | A clean checkout can reconstruct the test environment and CI exercises the relevant suite on every HT-branch push. |

## Frozen release gates

The branch must not be called HT-ready while any of these remain unresolved:

- exchange source units are incorrect or unverified;
- the setup validator rejects valid source-defined modes;
- a non-converged SCF can appear as successful;
- failed execution can lack persistent failure evidence;
- `reuse_if_identical` can reuse altered scientific inputs;
- symmetry reduction removes jobs outside a proven magnetic-reference domain;
- a clean-checkout test environment cannot be reconstructed;
- HT-R09 negative controls fail to detect intentionally corrupted cases.

## HT-R00 decision and evidence

HT-18 is preserved as historical evidence only. The readiness report now
states prominently that this branch is not ready for scientific HT
production, and this ledger freezes the audit findings for R01-R09.

No production code, parser, validator, CI workflow, or test was changed for
HT-R00. The remediation start commit is
`4981cbd952559df6f99fbc933e2091722e95f4db`, and the remediation start date is
`2026-08-31`.

## HT-R01 decision and evidence

The current [`fable_v4/source/exchange.f90`](https://raw.githubusercontent.com/rslmtoasa/rslmtoasa/fable_v4/source/exchange.f90)
audit establishes the source unit for the interaction writer paths. It writes
`jijso.out`/`jijfo.out`,
`dijso.out`/`dijfo.out`, and `aijso.out`/`aijfo.out` after multiplying their
calculated scalar, DMI, and full `aij` tensor values by
`1.0d3/(4.0d0*pi)`. The corresponding unsuffixed `jij.out`, `dij.out`, and
`aij.out` paths use the same scaling. These source values are mRy; the raw
numeric values are not rescaled by RS-utils.

RS-utils now exposes `Unit.MILLI_RY`/`MILLI_RYDBERG`, tags audited scalar,
DMI, and `aij` tensor outputs as `mRy`, and supports explicit Ry↔mRy
conversion without changing the parsed source quantity. Generic tensor
files without a source-unit declaration remain `unknown` rather than
receiving a guessed unit. Existing result schema version `v1` remains
backward-readable because each quantity already carries its unit; old records
declared as `Ry` are not silently migrated.

Regression coverage verifies unchanged MnSi source numbers, mRy JSON
round-tripping, the explicit 1000:1 conversion relation, source-backed tensor
units, and unknown units for an unverified tensor filename. The focused
HT-12, HT-11, and HT-14 suites and the full suite pass.

## Scope and dependency order

The existing worker-subprocess containment, transactional setup service,
typed core contracts, canonical `PairKey(i,j,R)` model, exchange job
ownership, parser framework, ASE structure boundary, TB model, and lightweight
campaign architecture remain in scope for focused repair. They are not to be
discarded unless a concrete defect requires a local change.

```text
HT-R00
  |
  +--> HT-R01
  +--> HT-R02
  +--> HT-R03
  +--> HT-R04
  +--> HT-R05
  +--> HT-R06
  +--> HT-R07
  +--> HT-R08
          \\        /
           \\      /
            HT-R09
```

The governing definition of a successful remediation is conservative:
supported domains must be scientifically correct and explicit, unsupported
or unverified domains must fail closed or be clearly marked, and legacy
valid workflows must remain available. Determinism, source units, identity,
execution/parse/convergence separation, tamper detection, conservative
symmetry, reproducible tests, CI coverage, and adversarial closeout evidence
are all required.

## HT-R02 source-derived validation contract

The current public [`fable_v4` lattice namelist declarations](https://github.com/rslmtoasa/rslmtoasa/blob/fable_v4/source/include_codes/namelists/lattice.f90),
[`lattice` implementation](https://github.com/rslmtoasa/rslmtoasa/blob/fable_v4/source/lattice_lifecycle.f90),
and [`charge` implementation](https://github.com/rslmtoasa/rslmtoasa/blob/fable_v4/source/charge.f90)
define these fields independently:

| Field | Source meaning used by validation |
| --- | --- |
| `nbulk_bulk` | Number of primitive-cell/reference bulk sites used by the mode-specific lattice setup. It participates in surface/interface derived counts; it is not a universal array length. |
| `ntot` | Number of sites in the primitive cell. `crd`, `izp`, `no`, `iu`, and `ib` are the primitive-cell tables and are sized by this count. |
| `nbas` | Number of shells/sites used by charge distribution. It can differ from `ntot`; the `charge` source explicitly allows `nbas` to differ from `ntot` for reduced cells and uses it for surface charge storage. |
| `nrec` | Number of sites/classes included in self-consistent recursion. `irec` is sized by this count. Bulk presets set it to the bulk count, while surface/interface construction derives it from the active layer count. |

The source also documents `ntype = nbulk + nrec`. For surface mode
(`calctype='S'`), `build_clusup` derives `ntype = nbulk_bulk + nlay` and
`nrec = nlay`. For layered/interface mode (`calctype='L'`), it derives
`ntype = 2*nbulk_bulk + nlay` and `nrec = nlay`. These relationships are
validated only when the corresponding mode fields are present. The validator
does not infer a missing mode, invent an upper bound for cluster indices, or
require `nbas`, `ntot`, `nbulk_bulk`, and `nrec` to be equal.

The public `fable_v4` surface and impurity inputs are source-realistic
examples of why labels cannot be tied to `nrec`: the surface input has seven
labels with six active layers, while the hcp impurity input has three labels
with two primitive/reference sites. No checked-in public unequal-count
`lattice.nml` was found, so regression coverage uses a minimal fixture built
from those source-declared surface semantics and records the provenance here.
The protected Fortran kernel remains responsible for cluster construction,
mode-specific geometry, and any bounds that depend on the generated cluster.

## HT-R02 decision and evidence

The validator no longer applies a universal `count == nrec` rule. It now
uses `ntot` for the primitive-cell `crd`, `izp`, `no`, `iu`, and `ib` tables,
`nrec` for `irec`, and leaves `nbas` and `nbulk_bulk` independent except for
the source-backed surface/interface derivations. Input labels are checked
against `ntype` and are no longer incorrectly required to equal `nrec`.
Surface (`S`) and layered/interface (`L`) diagnostics explicitly check their
source-derived `nlay` relationships when the mode is declared. Cluster index
upper bounds remain delegated to the protected kernel because the cluster is
constructed after preflight.

The supercell transform now scales `ntot` and only scales the other count
fields when the source explicitly used them as the primitive-site count;
independent mode-specific counts are retained. This preserves the ordinary
bulk workflow while avoiding a new count aliasing invariant.

Regression coverage includes the bcc baseline, a source-derived surface
fixture with `nbulk_bulk=1`, `ntot=1`, `nbas=51`, `nrec=6`, and seven labels, a
public fable_v4-style impurity label/count case, malformed recursion-table
and surface-layer cases, plus the existing supercell/setup-service and CLI
coverage. The required focused suites pass:

```text
tests/test_setup_validation.py       13 passed
tests/test_setup_service.py          11 passed
tests/test_legacy_cli.py              8 passed
tests/test_fortran_setup_runner.py   4 passed
full feasible suite                  104 passed, 10 warnings
```

## HT-R03 decision and evidence

Campaign status now retains independent execution, output/parser, and
scientific-convergence axes. `COMPLETED_VALID` is intentionally narrow: it
requires an atomically persisted `rsutil-execution-v1` record with zero exit,
parser-valid output, and explicit SCF convergence. Non-converged and
convergence-unknown results remain collectable but are not successful
scientific results. Missing or identity-mismatched execution evidence is
reported as unknown/legacy evidence, and a current non-zero execution remains
failed even when a stale `report.out` parses successfully.

The generated Slurm path writes a running record followed by an atomic final
record containing UTC timestamps, command, working directory, host, campaign
identity, available Slurm IDs, return code, signal information when observed,
and success/failure status. The status and collection JSON/CSV outputs carry
all three axes; collection's successful view contains only explicit converged
successes.

Regression coverage includes the required converged, non-converged, unknown,
malformed, non-zero/no-output, stale-output, missing-record, and wrong-identity
cases. The required focused suites pass:

```text
tests/test_ht17_campaign.py  12 passed
tests/test_ht16_campaign.py   4 passed
tests/test_ht11_results.py    6 passed
```

The authoritative evidence is the public `fable_v4` source links recorded in
the contract above; the protected numerical Fortran kernel was not modified.

## HT-R04 decision and evidence

Exchange manifest schema `rsutil-exchange-v2` now treats the hashed
`spec.pair_keys` list as the only authoritative pair identity. The retained
top-level `canonical_pair_list`, `pair_keys`, and `legacy_ijpair` fields are
validated compatibility views; both the job reader and result parser reject
tampered or inconsistent views. Unambiguous v1 manifests remain readable,
but no duplicated view is trusted without consistency validation.

Campaign run manifests now use the `rsutil-campaign-run-v2` contract and
record every setup-generated file in `campaign.materialized_inputs`, including
its SHA-256 and size. `reuse_if_identical` verifies the requested entry and
manifest identity, the materialized-input contract, every current file hash,
the campaign integrity digest, and unexpected `.nml` files. It rejects missing, modified, or unrecorded
scientific inputs with `OWNED_FILE_MISSING`, `INPUT_MODIFIED`, or
`UNEXPECTED_OWNED_FILE` evidence. Runtime files such as reports and logs are
outside this identity set and do not prevent valid reuse.

Regression coverage includes tampered top-level and authoritative pair data,
modified and missing generated campaign inputs, valid reuse, and runtime
output added after generation. The required focused suites pass:

```text
tests/test_ht10_exchange_jobs.py      7 passed
tests/test_ht12_exchange_results.py  10 passed
tests/test_ht16_campaign.py           7 passed
tests/test_ht17_campaign.py          12 passed
full feasible suite                 120 passed, 10 warnings
```

## HT-R05 decision and evidence

The exchange reducer now has an explicit `MagneticSymmetryContext`. The only
reduction domain claimed is collinear ferromagnetism: without SOC, exact
chemical periodic operations are used; with SOC, operations must also
preserve the declared global magnetization axis as an axial vector, including
the determinant sign for improper operations. The legacy `--soc` path remains
the `+z` global-axis convention.

AFM, ferrimagnetic, generic noncollinear, magnetic-space-group, time-reversal,
and unknown contexts are represented but resolve to identity-only reduction.
They cannot silently inherit chemical FM equivalence. The canonical
`RSStructure` site labels/types are carried into symmetry mapping, so same-
element sites with distinct RS identity cannot be merged. Missing symmetry
analyzer evidence also retains all neighbors. The job specification records
the magnetic reference class so generated-job identity includes the reduction
domain.

The current public [`fable_v4` control contract`](https://raw.githubusercontent.com/rslmtoasa/rslmtoasa/fable_v4/source/control.f90#L51-L66)
distinguishes collinear and noncollinear modes and SOC combinations, while
[`exchange.f90`](https://raw.githubusercontent.com/rslmtoasa/rslmtoasa/fable_v4/source/exchange.f90#L148-L166)
consumes the explicit pair list; no upstream magnetic-space-group reduction
contract was inferred.

Regression coverage includes FM periodic reduction, analyzer-absent fallback,
SOC axial filtering, distinct same-element RS labels, exact pair-orbit
identity, and no-reduction controls for AFM/ferrimagnetic/noncollinear and
other unproven contexts. The required focused suites pass:

```text
tests/test_ht07_pair_identity.py       6 passed
tests/test_ht08_exchange_service.py    6 passed
tests/test_ht09_exchange_symmetry.py    8 passed
tests/test_ht10_exchange_jobs.py        7 passed
full feasible suite                   126 passed, 12 warnings
```

## HT-R06 decision and evidence

The required current public `fable_v4` source audit was completed against
commit `faf71320c14afd680ce6973fb48fe4b6f8a25529`. The previous HT-13 claim
that no conductivity or damping source/writer existed was disproved.

The audit found a KPM/Kubo-style conductivity writer in
`source/conductivity.f90::calculate_conductivity_tensor`, with controls in
`&calculation post_processing='conductivity'`, `&control cond_type`,
`cond_calctype`, `cond_ll`, and `random_vec_num`. The source emits total,
per-type, and orbital real/imaginary files for charge, spin, orbital,
accumulation, and torque operator modes. It writes energy-relative-to-Fermi,
real, and imaginary numeric columns, but no units, headers, component map,
operator/channel metadata, broadening record, or frequency axis. The active
prefactor also omits the computed volume; the physical-constant alternatives
are comments. These outputs remain `UNSUPPORTED_AMBIGUOUS`, not a guessed
conductivity result.

The audit also found `source/exchange_dynamics.f90::calculate_gilbert_damping`
and its `alldampings.out` / `damping-energy.out` writers, gated through
`&calculation do_damping=.true.` and SOC/`nsp=2` conditions. Pair tensor order,
energy-relative-to-Fermi output, pair indexing, and the separate experimental
moment-of-inertia files were traced. The damping files do not emit units,
normalization, Green-function route/broadening, SOC/torque provenance, or run
identity; the inertia source explicitly calls its kernel raw and experimental.
They therefore remain `UNSUPPORTED_AMBIGUOUS` and are not exposed as a
`DampingResult`. Phenomenological LLG `alpha`, relaxation rates, and generic
imaginary response were explicitly kept distinct from ab-initio Gilbert
damping.

`docs/HT_TRANSPORT_DAMPING_FORMATS.md` now contains the source-grounded
inventory, exact output/control details, unsupported/absent classifications,
and a proposed upstream machine-readable contract. `rsutil.results.transport`
records the audited evidence and recognizes the actual fable_v4 filenames,
while continuing to return no guessed values. Focused regression coverage
checks that conductivity and damping filenames are detected but remain
fail-closed, and that registry evidence names the audited source files.

Required verification:

```text
tests/test_ht13_transport_results.py  17 passed
tests/test_ht11_results.py              6 passed
full feasible suite                    132 passed, 12 warnings
```

## HT-R07 decision and evidence

HT-R07 was completed without changing the protected Fortran kernel or
redesigning the TB/ASE workflows. Setup provenance now derives reproducibility
from active stochastic operations: deterministic unseeded setup is
reproducible, while lattice/basis distortion and SQS operations identify their
seeds and warn when an active operation is unseeded. Setup manifests receive
operation start/finish timestamps from the service lifecycle, and those
timestamps remain outside deterministic scientific identity hashes.

The supported ASE boundary remains `Atoms <-> RSStructure`. The full
ASE-to-setup-source helper is now explicitly experimental and requires
`experimental=True`; its generated LMTO policy is documented as placeholder
and reviewed rather than silently advertised as HT-ready. Its derived
Wigner-Seitz radius uses `(3 V_atom / (4 pi))^(1/3)`.

`Workspace.child()` now resolves the owned root and candidate path before
containment checking, rejecting `..`, absolute, and symlink escapes. Exchange
neighbor reduction accepts an explicit `index_base` and normalizes once to
canonical zero-based site indices; no integer-domain heuristic remains. The
unreachable post-service setup implementation was removed from `cli.py` after
the existing service-path CLI regressions passed.

Legacy `rs2txt.dat` conversion now solves its Cartesian displacement equation
directly in the supplied lattice basis, preserving exact long-range periodic
identities instead of limiting the search to translations in `[-5, 5]`.
The legacy path remains fail-closed for singular, malformed, or ambiguous
geometry. Relevant setup, ASE, path, indexing, CLI, and TB documentation was
updated.

Required verification:

```text
HT-R07 focused suites                 55 passed
full feasible suite                  136 passed, 12 warnings
```

## HT-R08 decision and evidence

HT-R08 separates PEP 517 build requirements from package runtime and optional
feature dependencies. The build system now contains only setuptools, wheel,
CMake, f90wrap, and NumPy needed to build the extension. Core runtime
requirements remain NumPy, f90nml, and pymatgen; ASE, band/DOS dependencies,
and PyVista are exposed through named extras. The `test` extra includes every
dependency imported by the complete checked-in suite, including ASE and the
band-analysis stack. The package remains at version `0.10.0`; no release is
being cut by this task.

`docs/HT_TESTING.md` records the supported Python range and exact Linux
clean-checkout commands, including the gfortran/CMake prerequisites, editable
PEP 517 install, wheel build, and `python -m pytest -q` verification command.
The CLI visualization import is deferred until that feature is requested, so
the core package and setup/exchange CLI do not make PyVista mandatory.

The pip CI workflow now installs `.[test]`, installs the Linux Fortran
toolchain, builds a wheel, and runs the same `python -m pytest -q` suite on
every push and pull request. The cibuildwheel workflow is no longer excluded
from every branch; it runs for `main` and `HT` pushes and pull requests, with
wheel smoke tests enabled on the supported Linux toolchain. The native CMake
workflow remains the cross-platform build check. Cibuildwheel's Linux image
installs the required Fortran/CMake toolchain before building.

Local verification from the existing environment:

```text
python -m pytest --collect-only -q       137 tests collected
python -m pytest -q                      137 passed, 12 warnings
python -m build --wheel                 passed with isolated build dependencies
temporary venv wheel import              passed outside the source checkout
```

The clean-environment commands are documented and wired into CI. The built
wheel includes the compiled extension and imports successfully outside the
source checkout; the temporary wheel check used the existing installed runtime
libraries and did not claim a second dependency installation.

## HT-R09 decision and evidence

HT-R09 was completed as a falsification attempt against the R00-R08 release
gates. No protected Fortran numerical source, algorithm, constant, tolerance,
or formula was changed. The current readiness decision and support matrix are
in [`HT_READINESS_REPORT.md`](HT_READINESS_REPORT.md); the historical HT-18
evidence remains preserved there under its own heading.

### Reproducible full-suite environment

The closeout was run from the R08 environment using the repository's isolated
PEP 517 wheel-build path and the documented test extra. The observed local
environment was:

```text
Python 3.11.16       /Users/andersb/envs/p311/bin/python
platform             macOS-26.5.2-arm64-arm-64bit (arm64)
gfortran             GNU Fortran (Homebrew GCC 16.2.0)
cmake                4.4.3
install/build path   python -m build --wheel; temporary venv wheel import;
                     pip install -e '.[test]' for the checkout test run
```

`python -m pytest -q` passed with **137 passed, 12 warnings** in 57.51 s.
The warnings are the known spglib deprecation warnings listed by pytest. The
focused adversarial/legacy suites passed with **86 passed, 12 warnings** in
49.50 s, and the legacy CLI plus setup compatibility subset passed with
**22 passed** in 37.73 s. No remote CI service was available from this
environment; CI workflow configuration was inspected but no remote run is
claimed.

### Negative-control matrix

The focused suite and direct parser control covered the required failure
boundaries:

| Control | Evidence | Result |
| --- | --- | --- |
| Ry/mRy mismatch or mixed declaration | `tests/test_ht12_exchange_results.py`; direct `jijso.out` mixed `Ry mRy` parser control | conflicting units rejected; conversion is explicit and `1 mRy = 0.001 Ry` |
| Valid unequal counts / inconsistent counts | `tests/test_setup_validation.py::test_source_derived_surface_counts_can_differ`; `...::test_recursion_array_count_is_still_rejected`; `...::test_inconsistent_ntype_is_coded_as_a_semantic_error` | pass / rejected |
| Exit 0 + `Not converged` | `tests/test_ht17_campaign.py::test_nonconverged_report_is_not_a_scientific_success` | non-converged |
| Nonzero exit with stale good report | `tests/test_ht17_campaign.py::test_nonzero_exit_cannot_be_hidden_by_stale_parseable_report` | execution failure |
| Malformed report | `tests/test_ht17_campaign.py::test_zero_exit_with_malformed_report_is_parser_failure` | parser invalid |
| Unknown convergence | `tests/test_ht17_campaign.py::test_zero_exit_with_unknown_convergence_remains_explicit` | unknown, never success |
| Reuse unchanged / mutated / missing materialized input | `tests/test_ht16_campaign.py` reuse parametrization and missing-input tests | unchanged reused; mutation/missing rejected |
| Exchange duplicate, stale, or unauthenticated pair mapping | `tests/test_ht10_exchange_jobs.py` and `tests/test_ht12_exchange_results.py` tamper tests | rejected without replacement |
| AFM/noncollinear symmetry; RS identity | `tests/test_ht08_exchange_service.py`, `tests/test_ht09_exchange_symmetry.py` | no unproven reduction; site/type identity preserved |
| Worker `STOP`/failure isolation and partial publication | `tests/test_fortran_setup_runner.py` | parent survives; failed destination is not published |
| Deterministic vs active stochastic unseeded setup | `tests/test_setup_validation.py` deterministic/stochastic manifest tests | deterministic reproducible; active stochastic marked nonreproducible |
| Path escape and ambiguous indexing | `tests/test_core_contracts.py`; `tests/test_ht09_exchange_symmetry.py::test_reference_indices_require_an_explicit_index_base` | rejected |

### Legacy compatibility

`tests/test_legacy_cli.py` and the setup-validation compatibility cases passed
for `rsutil -S`, `rsutil -J`, simple bulk, multisite/multispecies, seeded
supercell, exchange jobs, and a source-derived valid unequal-count surface
mode. Compatibility is semantic/workflow compatibility, not byte identity.

### Independent pilot evidence

A fresh controlled pilot was generated in two independent temporary roots.
It contained bcc Fe, FeN multisite/multispecies, FePt SOC, seeded bcc
`2x1x1` supercell/distortion, and a valid source-derived unequal-count
surface mode. The real Fortran setup runner generated the first four cases;
the surface fixture used a controlled setup runner because its purpose was to
exercise the validated count semantics without asserting a scientific surface
potential. The calculation stage used a deliberately controlled executable,
not a scientific RS-LMTO run.

The `lattice_new.nml` and `input_new.nml` SHA-256 hashes matched between both
roots for every entry. The pipeline was exercised as generate -> render
`submit_array.sbatch` -> execute controlled fixtures -> status -> collect.
The controlled outcomes in each root were:

```text
returns             [0, 0, 0, 9, 0]
statuses             COMPLETED_VALID, COMPLETED_NOT_CONVERGED,
                     COMPLETED_INVALID, FAILED, COMPLETED_UNKNOWN
scientific states    CONVERGED, NOT_CONVERGED, UNKNOWN,
                     NOT_APPLICABLE, UNKNOWN
scheduler state      not_queried
collected success    bcc-fe only
```

The converged, non-converged, malformed, nonzero-exit, and unknown fixtures
were intentional controls. They are evidence for state separation and
collection behavior, not evidence of scientific convergence. No scheduler
execution was claimed.

### R09 gate decision

R00-R08 are checked and their evidence sections were reread before this
closeout. The full suite, focused negative controls, legacy compatibility
subset, and independent deterministic pilot hash comparison passed. Therefore
the foundation is **ready for pilot use** within the support matrix in
`HT_READINESS_REPORT.md`; it is not a blanket production-readiness claim.

## Final closeout HT-R10…HT-R14

The earlier R09 conclusion is **superseded/reopened by the post-R09
self-review**. Its historical evidence remains above and in
[`HT_READINESS_REPORT.md`](HT_READINESS_REPORT.md); it is not current
readiness evidence. The final closeout evidence below supersedes that
historical decision within the stated support matrix.

- [x] **HT-R10** — Interaction scientific contract closeout
- [x] **HT-R11** — Authoritative regeneration-based reuse
- [x] **HT-R12** — Magnetic-reference-aware symmetry contract
- [x] **HT-R13** — Installed-wheel clean-environment verification
- [x] **HT-R14** — Final falsification closeout

### HT-R10 decision and evidence

The current public `fable_v4` `source/exchange.f90` writer contract was
audited directly. `jij`, `dij`, and `aij` outputs are scaled by
`1.0d3/(4.0d0*pi)`, so their source unit is mRy. The writer emits `aij` as
the anisotropic component; the complete tensor is constructed separately as
`J*I + [[0,Dz,-Dy],[-Dz,0,Dx],[Dy,-Dx,0]] + A`. Whole-array Fortran output is
column-major, and the parser now reconstructs that order explicitly.

Known `fable_v4` filenames now reject contradictory explicit unit headers with
`SOURCE_UNIT_MISMATCH`; absent headers retain the source-defined mRy unit.
Generic formats remain explicit-header/unknown-unit paths. `ExchangePair`
geometry carries a length unit through result and `R_ij` JSON/text
serialization, with legacy unresolved data labeled `unknown` rather than
assumed to be Angstrom.

The raw tensor model is `anisotropic_exchange_tensor`; the former
`exchange_tensor`/`tensor` attributes remain deprecated read aliases. The
derived `full_exchange_tensor()` operation requires J, D, and A with known,
compatible energy units and otherwise reports unavailable.

Focused evidence:

```text
python -m pytest -q tests/test_ht12_exchange_results.py tests/test_ht07_pair_identity.py
23 passed
python -m pytest -q tests/test_ht12_exchange_results.py tests/test_ht07_pair_identity.py tests/test_ht14_tight_binding.py
28 passed
python -m pytest -q
144 passed, 12 warnings
```

### HT-R11 decision and evidence

HT-R11 is closed for campaign and exchange reuse. Campaign reuse now
regenerates each setup in an isolated candidate directory from the current
source and `CampaignEntry` setup specification, then compares the candidate's
immutable output hashes and stable setup identity with the existing run.
Exchange reuse already stages a fresh candidate; it now compares destination
files against that candidate's hashes rather than hashes read from the old
destination manifest. Both paths reject unseeded stochastic generation and
allow runtime artifacts outside the immutable generated-file contract. The
campaign path rejects unseeded stochastic generation conservatively.

The historical manifest integrity hash remains useful as tamper evidence, but
is not a root of trust: a manifest can self-consistently record a coordinated
file edit. The fresh-regeneration contract is recorded as
`fresh_regeneration_v1` in campaign manifests and is documented in
`HT_CAMPAIGNS.md` and `HT_EXCHANGE_JOBS.md`.

Focused R11 evidence includes negative controls that modify both a generated
scientific file and its historical hash record; both are rejected against
freshly generated expectations. Positive controls cover campaign and exchange
runtime outputs and seeded stochastic campaign reuse:

```text
python -m pytest -q tests/test_ht10_exchange_jobs.py tests/test_ht16_campaign.py
21 passed, 1 warning
```

### HT-R12 decision and evidence

HT-R12 is closed for source-derived magnetic context. The current public
`fable_v4` `control` source documents `nsp=1` as collinear scalar-relativistic,
`nsp=2` as collinear fully-relativistic, `nsp=3` as noncollinear
scalar-relativistic, and `nsp=4` as noncollinear fully-relativistic. Its
`is_collinear`, `is_noncollinear`, and `has_soc` helpers implement exactly
that mapping. These values describe the spin representation and SOC mode;
they do not establish FM, AFM, or ferrimagnetic order.

`ReferenceCalculation` now retains `control.nsp` and exposes the centralized
source-derived `nsp`, collinearity, noncollinearity, and SOC properties.
Missing and unsupported values remain explicit unknowns. Exchange planning
derives SOC from the reference when present and rejects an explicit conflicting
SOC argument rather than silently overriding the source.

Symmetry planning defaults to the `unknown` magnetic-reference class. Unknown
order, AFM, ferrimagnetic, and noncollinear references retain the complete
candidate pair set. Only an explicit `collinear_fm` assertion can enable the
proven reduction domain; SOC additionally requires an explicit global
magnetization axis. The CLI remains compatible with `rsutil -J`, and now
offers `--magnetic-reference collinear_fm --magnetic-axis X Y Z` for that
restricted opt-in. Job manifests record source `nsp`, derived spin/SOC state,
magnetic class, axis, reduction status, and the conservative-reduction reason.

Negative controls cover all four `nsp` modes, missing `nsp`, conflicting SOC,
noncollinear references, absent SOC axes, and accidental FM inference.
Positive controls cover explicit collinear-FM reduction with and without SOC,
including the restricted axial-axis case. Same-element RS labels remain part
of site identity and the exact periodic pair-key tests remain unchanged.

Focused R12 evidence:

```text
python -m pytest -q tests/test_ht12_magnetic_reference.py tests/test_ht08_exchange_service.py tests/test_ht09_exchange_symmetry.py tests/test_ht10_exchange_jobs.py tests/test_legacy_cli.py
45 passed
```

### HT-R13 decision and evidence

HT-R13 is closed for installed-wheel clean-environment verification. CI now
has separate source/editable and artifact lanes. The artifact lane builds
exactly one wheel, installs it with the declared `test` extra into an
external venv, checks dependencies, proves both the Python package and
compiled extension resolve under that venv's site-packages, and runs the
complete suite from an external copy of `tests` and `examples` with
`PYTHONPATH` removed. The workflow runs on both pushes and pull requests.

The artifact run also restored the legacy CLI's documented default recursion
depth (`lld=21` for the default method); the existing legacy test caught the
incompatibility when the installed wheel, rather than an unrelated global
`rsutil` executable, was first placed on `PATH`.

Focused R13 evidence from a clean local wheel build:

```text
platform: macOS-26.5.2-arm64-arm-64bit (arm64)
Python: 3.11.16
gfortran: GNU Fortran (Homebrew GCC 16.2.0)
cmake: 4.4.3
wheel: rsutil-0.10.0-cp311-cp311-macosx_26_0_arm64.whl
rsutil.__file__: .../site-packages/rsutil/__init__.py
native extension: .../site-packages/rsutil/rsutil.so
site-packages assertion: PASS
pip check: No broken requirements found
pytest: 162 passed, 15 warnings
```

### HT-R14 decision and evidence

HT-R14 is closed. The final falsification pass was evaluated against code
commit `8840afa6b0ac0ef707e2a414d8d205112faff866`, with the documentation
closeout committed separately. The source checkout passed 162 tests with 15
known spglib deprecation/user warnings. The clean external wheel environment
also passed 162 tests with the same 15 warnings; `pip check` was clean and
both the Python package and native extension resolved from site-packages.

The mandatory negative controls were exercised by the focused suites and the
full artifact suite: contradictory known-format units, raw-versus-derived
tensor semantics and Fortran ordering, explicit pair geometry units, valid
and malformed unequal-count setup relationships, every execution/convergence
state, coordinated campaign/exchange tampering, pair-manifest tampering,
all four `nsp` modes and conservative magnetic reduction, worker failure and
transactional publication, seeded/unseeded reproducibility, path/index
boundaries, and installed-artifact import provenance. No repaired contract
collapsed into a misleading success.

The deterministic `examples/campaigns/pilot.json` campaign was generated in
two independent temporary roots. Entry IDs, generated scientific input
hashes, and scientific manifest fields matched for all three entries,
including the seeded SQS case. The existing HT-17 controlled execution
fixtures cover converged, non-converged, unknown, malformed, nonzero-exit,
stale-report, and signaled outcomes. Legacy `rsutil -S` and `rsutil -J`
regressions passed semantically through the installed-wheel suite.

No remote CI result was available from this checkout: `origin/HT` remains at
`54a279b6e00b5a51009626d94a16d1964149fc45`, before the closeout commits. This
is recorded as unavailable evidence, not as a remote pass claim. The support
boundary and remaining unsupported domains are recorded in
`docs/HT_READINESS_REPORT.md`.
