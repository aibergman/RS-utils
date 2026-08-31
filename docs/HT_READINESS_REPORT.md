# RS-utils HT readiness report — HT-R15…R17 final shape-up

> **READY FOR PILOT HT USE WITHIN THE DOCUMENTED SUPPORT MATRIX — HT-R17
> COMPLETE.**

The R15/R16 identity, ownership, and magnetic-context shape-up is complete.
The installed-wheel closeout passed 177 tests with 16 known spglib deprecation
warnings. The wheel was `rsutil-0.10.0-cp311-cp311-macosx_26_0_arm64.whl`;
both the Python package and compiled extension resolved from the external
environment's site-packages with `PYTHONPATH` removed. The source checkout
suite passed the same 177 tests and warning count.
>
> The earlier R09 conclusion is superseded/reopened by the post-R09
> self-review. Its historical evidence is retained below, but it must not be
> read as universal RS-LMTO feature coverage. The protected Fortran
> numerical kernel was not modified.

## Historical R09 closeout (superseded)

### Reproducible environment and test result

The closeout ran from the R08 environment:

```text
Python 3.11.16       /Users/andersb/envs/p311/bin/python
platform             macOS-26.5.2-arm64-arm-64bit (arm64)
gfortran             GNU Fortran (Homebrew GCC 16.2.0)
cmake                4.4.3
build/install        python -m build --wheel; temporary venv wheel import;
                     pip install -e '.[test]' for the checkout test run
```

The required full command, `python -m pytest -q`, passed: **137 passed, 12
warnings** in 57.51 s. The warnings are known spglib deprecation warnings.
The focused adversarial suites passed: **86 passed, 12 warnings** in 49.50 s;
the legacy CLI/setup subset passed: **22 passed** in 37.73 s. No remote CI
service was available, so no remote CI run is claimed; CI workflow coverage is
documented by R08.

### Negative controls

| Contract boundary | Result and evidence |
| --- | --- |
| Units: wrong/mixed Ry vs mRy | Explicit conversion is required; `1 mRy -> 0.001 Ry`. A mixed `Ry mRy` `jijso.out` declaration was rejected. See `tests/test_ht12_exchange_results.py`. |
| Unequal setup counts | Valid source-derived surface counts passed; missing recursion entries and inconsistent `ntype` were rejected. See `tests/test_setup_validation.py`. |
| Execution/state separation | Exit 0 + `Not converged` => `COMPLETED_NOT_CONVERGED`; nonzero exit, including stale good output => `FAILED`; malformed => `COMPLETED_INVALID`; unknown convergence never succeeds. See `tests/test_ht17_campaign.py`. |
| Reuse and identity | Unchanged materialized input reused; mutation or missing input rejected. See `tests/test_ht16_campaign.py`. |
| Exchange mapping | Source-faithful `izp` site-to-type mapping, repeated-type geometry mapping, and duplicate/stale/tampered pair views are covered; exchange output type columns are never treated as site numbers. See `tests/test_ht15_site_types.py`, `tests/test_ht10_exchange_jobs.py`, and `tests/test_ht12_exchange_results.py`. |
| Magnetic symmetry and RS identity | AFM/noncollinear cases do not receive unproven reduction; site/type and distinct RS labels remain identity-bearing. See `tests/test_ht08_exchange_service.py` and `tests/test_ht09_exchange_symmetry.py`. |
| Worker failure and publication | Worker failure/`STOP` is isolated from the parent and partial destinations are not published. See `tests/test_fortran_setup_runner.py`. |
| Randomness | Unseeded deterministic setup is reproducible; active unseeded stochastic setup is marked nonreproducible. See `tests/test_setup_validation.py`. |
| Paths and indexing | Absolute/`..`/symlink escapes and implicit index-base ambiguity are rejected. See `tests/test_core_contracts.py` and `tests/test_ht09_exchange_symmetry.py`. |

### Legacy compatibility

The legacy CLI/setup subset passed for `rsutil -S`, `rsutil -J`, simple bulk,
multisite/multispecies, seeded supercell, exchange jobs, and a valid
source-derived unequal-count surface mode. Compatibility is semantic and
workflow-oriented, not byte-for-byte identity.

### Independent controlled pilot

The pilot was regenerated from scratch in two independent temporary roots. It
contained bcc Fe, FeN multisite/multispecies, FePt SOC, seeded bcc `2x1x1`
supercell/distortion, and a valid source-derived unequal-count surface mode.
The real Fortran setup runner generated the first four cases. The surface
fixture used a controlled setup runner to exercise count semantics without
claiming a scientific surface-potential result. The calculation executable was
deliberately controlled and is not a scientific RS-LMTO run.

| Entry | `lattice_new.nml` SHA-256 | `input_new.nml` SHA-256 |
| --- | --- | --- |
| bcc-fe | `cc100ed40943999ef6c71c82f80af1c0d64bf7bc63bfc4e396e0da31db9fdb4d` | `482b23fcdde989a98ebb7a538f9b65b87166106f9d4844aebe362a99510cb41b` |
| fen-multisite | `f38d20539652cf184c4932c86a3cdf2542f9c280ff2ab4d7221ef36345e80039` | `24e6536c6cab38eb99ff48730d65833b627667c57898070e1950dc66f4e183c0` |
| fept-soc | `59c618d67a9c54c9060d001c8c78fc9d0bb97d0b7b91a5bde08f93a7ff0a36cd` | `038e56451fe05600347d8d99927da87c57b585124d119692083210503f9baf9a` |
| bcc-seeded-supercell | `8e57ac5fe9efed5ff0e89fb8f9225cc4b48839a56d1e4b5910a442b1783bd6e4` | `4013b633cce85fa7a890b6016f256c1267a3f09befc864c86bb3a723d37efd3f` |
| surface-unequal-counts | `186c3026e9dd288d554919b7d159abf809b41eee5fc9e71492756ab5828388d5` | `ed6e1819c39e9ffabf95ee74ffac73dfc67a2389a0093ca3f8bcf6f7a196148e` |

Both roots matched byte hashes for every listed scientific input. The
pipeline was exercised as `generate -> render -> fake/controlled execute ->
status -> collect`; the rendered scheduler state remained `not_queried`.
Controlled returns were `[0, 0, 0, 9, 0]`, producing respectively
`COMPLETED_VALID`, `COMPLETED_NOT_CONVERGED`, `COMPLETED_INVALID`, `FAILED`,
and `COMPLETED_UNKNOWN`. Only bcc Fe was collected as successful. These are
state-machine controls, not scientific convergence evidence.

## Current support matrix

Statuses are deliberately limited to `READY`, `RESTRICTED`, `EXPERIMENTAL`,
and `UNSUPPORTED`.

| Capability | Status | Evidence and boundary |
| --- | --- | --- |
| Ordinary bulk setup | READY | Real Fortran worker and bcc/FeN/FePt fixtures pass; isolated and transactional; `tests/test_fortran_setup_runner.py`, `tests/test_setup_validation.py`. |
| HCP setup | RESTRICTED | Generic lattice validation is covered, but no checked-in HCP end-to-end worker fixture is part of this closeout; `tests/test_setup_validation.py`. |
| Supercell setup | READY | Seeded supercell generation is deterministic and non-destructive; `tests/test_legacy_cli.py`, `tests/test_setup_validation.py`. |
| Surface setup | RESTRICTED | Source-derived unequal-count semantics validate; surface potential generation is not claimed; `tests/test_setup_validation.py`. |
| Impurity/defect setup | RESTRICTED | Source-derived impurity label/count semantics validate; no broad defect workflow claim; `tests/test_setup_validation.py`. |
| Stochastic transformations/SQS | READY | Seeded transformations reproduce; active unseeded stochastic generation is rejected; `tests/test_setup_validation.py`, `tests/test_legacy_cli.py`, `tests/test_ht16_campaign.py`. |
| Scalar exchange J | READY | Canonical pair generation, manifests, and real MnSi result parsing preserve source mRy; `tests/test_ht12_exchange_results.py`, `tests/test_ht07_pair_identity.py`. |
| DMI vector D | READY | Real MnSi `dijso.out` parses with source mRy provenance; `tests/test_ht12_exchange_results.py`. |
| Raw anisotropic A tensor | EXPERIMENTAL | Source-faithful meaning and ordering are tested, but unverified external tensor sources retain unknown units; `tests/test_ht12_exchange_results.py`. |
| Derived full exchange tensor | EXPERIMENTAL | Exact `J I + D(antisym) + A` derivation and missing-component failure are tested; broader production tensor coverage remains bounded; `tests/test_ht12_exchange_results.py`. |
| Pair geometry | READY | Distance and displacement retain explicit length units/conventions through JSON/text export; `tests/test_ht12_exchange_results.py`, `tests/test_legacy_cli.py`. |
| No-SOC FM symmetry | READY | Exact reduction is enabled only by explicit collinear-FM context; `tests/test_ht12_magnetic_reference.py`, `tests/test_ht09_exchange_symmetry.py`. |
| SOC FM symmetry | RESTRICTED | Source `nsp` controls SOC; reduction additionally requires an explicit collinear-FM reference and axial magnetization axis; `tests/test_ht12_magnetic_reference.py`, `tests/test_ht09_exchange_symmetry.py`. |
| AFM/ferrimagnetic symmetry | RESTRICTED | Conservative identity/full-pair behavior is covered; no general magnetic reduction is claimed; `tests/test_ht08_exchange_service.py`, `tests/test_ht09_exchange_symmetry.py`. |
| Generic noncollinear symmetry | UNSUPPORTED | No general equivalence theorem is asserted; full inventory and domain review are required. |
| SCF parsing/convergence | READY | Converged, non-converged, malformed, unknown, and execution-failure states are distinct and tested. |
| Conductivity | UNSUPPORTED | Source audit found insufficient machine-readable unit/operator context; parser fails closed. |
| Damping | UNSUPPORTED | Source audit found insufficient validated format/context; no scientific readiness claim. |
| Legacy interaction output formats | RESTRICTED | Audited scalar/DMI current formats are ready; legacy or unverified tensor variants require explicit review; `tests/test_ht12_exchange_results.py`. |
| Versioned current interaction formats | READY | Versioned result schemas preserve source units, provenance, and tensor meaning for the supported scalar/DMI paths; `tests/test_ht12_exchange_results.py`. |
| Versioned TB | READY | Versioned metadata, index base, units, provenance, and export tests pass. |
| Legacy TB | RESTRICTED | Exact conversion/fail-closed geometry checks pass; legacy source formats require explicit review. |
| ASE structure interop | READY | `Atoms <-> RSStructure` conversion and identity/unit contracts pass. |
| ASE automatic setup | EXPERIMENTAL | Requires explicit opt-in and is documented as placeholder policy, not production setup generation. |
| Campaign generation | READY | Deterministic identity, collision, seeded reproducibility, and independent-root pilot checks pass; `tests/test_ht16_campaign.py`. |
| Authoritative campaign/exchange reuse | READY | Fresh regeneration rejects scientific-file plus historical-manifest coordinated tampering; `tests/test_ht10_exchange_jobs.py`, `tests/test_ht16_campaign.py`. |
| Slurm rendering | READY | Shell-safe array script/run-list rendering is tested locally. |
| Slurm execution evidence | RESTRICTED | Controlled local execution evidence passes; no scheduler or remote cluster execution was available. |
| Convergence-aware status/collection | READY | Filesystem-only status, parse failures, sidecars, hashes, continuation, stale-report, and execution-state behavior pass; `tests/test_ht17_campaign.py`. |
| Collection | READY | Controlled result collection and sidecar/provenance handling are covered; `tests/test_ht17_campaign.py`. |
| AiiDA integration | UNSUPPORTED | Outside the validated R00-R09 foundation boundary. |
| UppASD workflow integration | UNSUPPORTED | Outside the validated R00-R09 foundation boundary. |
| Database backend | UNSUPPORTED | No database backend is implemented or claimed in this foundation. |
| ML training | UNSUPPORTED | No model-training pipeline is implemented or claimed in this foundation. |

The R09 gate was recorded as satisfied for pilot use, but that conclusion was
superseded/reopened by the post-R09 self-review. R10-R14 are now closed for
pilot use within this matrix. The historical HT-18 claims below are retained
for traceability only.

## Current HT-R10 closeout

HT-R10 is closed for the interaction scientific contract. The parser now
uses the audited `fable_v4` mRy writer contract for known `jij*`, `dij*`, and
`aij*` files and rejects contradictory explicit unit headers with
`SOURCE_UNIT_MISMATCH`. Raw `aij` data is represented as an anisotropic
tensor in Fortran column-major order; `full_exchange_tensor()` derives the
complete `J*I + skew(D) + A` tensor only when all three known, compatible
quantities are present. Pair displacement and distance exports carry an
explicit length unit.

The focused R10/parser compatibility suites pass:

```text
python -m pytest -q tests/test_ht12_exchange_results.py tests/test_ht07_pair_identity.py tests/test_ht14_tight_binding.py
28 passed
python -m pytest -q
144 passed, 12 warnings
```

This section records the R10 evidence; the later R12-R14 sections record the
remaining closeout work and final decision.

## Current HT-R11 closeout

HT-R11 is closed for authoritative reuse. Campaign reuse regenerates the
current setup in an isolated candidate directory and compares immutable output
bytes against that fresh candidate. Exchange reuse compares existing immutable
files against the freshly staged job, while allowing runtime `*.out`, `*.log`,
and `rsutil_execution.json` artifacts. Historical manifest hashes are not a
root of trust, because a coordinated edit can make a mutable manifest agree
with altered bytes; the fresh-regeneration contract is explicit in the
campaign and exchange documentation.

The focused adversarial campaign/exchange suite passed:

```text
python -m pytest -q tests/test_ht10_exchange_jobs.py tests/test_ht16_campaign.py
21 passed, 1 warning
```

The final R14 decision is recorded below.

## Current HT-R12 closeout

HT-R12 is closed for the magnetic-reference-aware symmetry contract. The
public `fable_v4` `control` source establishes `nsp=1/2` as collinear and
`nsp=3/4` as noncollinear, with SOC in `nsp=2/4`. RS-utils now retains and
propagates that source value; it does not infer FM from collinearity. Unknown,
AFM, ferrimagnetic, and noncollinear references use identity-only reduction.
An explicit collinear-FM assertion is required for reduction, and SOC
reduction additionally requires an explicit global magnetization axis.

The focused R12 suite passed:

```text
python -m pytest -q tests/test_ht12_magnetic_reference.py tests/test_ht08_exchange_service.py tests/test_ht09_exchange_symmetry.py tests/test_ht10_exchange_jobs.py tests/test_legacy_cli.py
45 passed
```

R14 is addressed by the final closeout below.

## Current HT-R13 closeout

HT-R13 is closed for installed-wheel clean-environment verification. CI now
separates the source/editable test lane from the built-wheel artifact lane.
The artifact lane builds exactly one wheel, installs it with the declared test
dependencies into an external venv, removes `PYTHONPATH`, asserts that both
`rsutil.__file__` and the compiled extension resolve under site-packages, and
runs the full suite from an external copy of the test fixtures.

The clean local artifact verification passed with Python 3.11.16,
gfortran 16.2.0, and CMake 4.4.3:

```text
rsutil.__file__: .../site-packages/rsutil/__init__.py
native extension: .../site-packages/rsutil/rsutil.so
site-packages assertion: PASS
pip check: No broken requirements found
pytest: 162 passed, 15 warnings
```

## Current HT-R14 closeout

HT-R14 is closed for pilot use within the support matrix above. The final
falsification pass was evaluated against code commit
`8840afa6b0ac0ef707e2a414d8d205112faff866`; the closeout documentation is
committed separately.

### Installed-artifact and full-suite evidence

The exact wheel was
`rsutil-0.10.0-cp311-cp311-macosx_26_0_arm64.whl`. It was installed with
the `test` extra into an external Python 3.11.16 virtual environment, with
`PYTHONPATH` removed and tests/examples copied outside the checkout. Both
`rsutil.__file__` and `rsutil.rsutil` resolved from that environment's
`site-packages`; `pip check` passed. The external full suite passed **164
tests, 15 warnings**. The source-checkout suite independently passed the same
164 tests and 15 warnings.

The follow-up artifact paths observed in that run were:

```text
rsutil.__file__: /private/tmp/rsutil-r14-venv.BVRc7L/lib/python3.11/site-packages/rsutil/__init__.py
native extension: /private/tmp/rsutil-r14-venv.BVRc7L/lib/python3.11/site-packages/rsutil/rsutil.so
site-packages assertion: PASS
```

The warnings are the known spglib deprecation/user warnings recorded by the
focused symmetry tests; no test was skipped or suppressed for this closeout.
The local toolchain was Python 3.11.16, gfortran 16.2.0, and CMake 4.4.3 on
macOS arm64.

### Negative-control summary

| Control | Result | Evidence |
| --- | --- | --- |
| C1 interaction units | Contradictory known-format `Ry` headers reject; mRy remains source-preserved; conversion is factor 1000. | `tests/test_ht12_exchange_results.py` |
| C2 tensor semantics | Raw `aij` remains distinct from derived full tensor; distinct nine-value ordering and missing-component failures are covered. | `tests/test_ht12_exchange_results.py` |
| C3 pair geometry | Distance/displacement exports retain explicit length units. | `tests/test_ht12_exchange_results.py`, `tests/test_legacy_cli.py` |
| C4 setup semantics | Valid unequal-count surface/impurity relationships pass; malformed counts fail. | `tests/test_setup_validation.py` |
| C5 execution/convergence | Zero/converged, zero/non-converged, zero/unknown, malformed, nonzero, stale-report, and signaled cases remain distinct. | `tests/test_ht17_campaign.py` |
| C6 authoritative reuse | Scientific-file edits, coordinated with historical manifest/hash edits, are rejected; runtime outputs and seeded stochastic reuse pass; unseeded stochastic reuse is rejected. | `tests/test_ht10_exchange_jobs.py`, `tests/test_ht16_campaign.py` |
| C7 pair manifest | Canonical/compatibility pair-list tampering is rejected. | `tests/test_ht10_exchange_jobs.py`, `tests/test_ht12_exchange_results.py` |
| C8 magnetic reference | `nsp=1..4`, missing `nsp`, conflicting SOC, noncollinearity, absent axes, and accidental FM inference are covered conservatively. | `tests/test_ht12_magnetic_reference.py`, `tests/test_ht09_exchange_symmetry.py` |
| C9 Fortran containment | Worker fatal exit leaves the parent alive; failed candidates are not published. | `tests/test_fortran_setup_runner.py`, `tests/test_ht10_exchange_jobs.py` |
| C10 provenance/randomness | Deterministic setup has reproducible identity without a meaningless seed; active unseeded stochastic setup is non-reproducible. | `tests/test_setup_validation.py` |
| C11 paths/indexes | Unsafe paths and ambiguous index bases are rejected. | `tests/test_core_contracts.py`, `tests/test_ht09_exchange_symmetry.py` |
| C12 installed artifact | Package and compiled extension resolve from external site-packages with no source shadowing. | R13 wheel lane; external full-suite run above |

### Pilot and legacy compatibility

`examples/campaigns/pilot.json` was generated twice into independent
temporary roots. All three entry IDs, `lattice_new.nml` hashes,
`input_new.nml` hashes, and scientific manifest fields matched, including the
seeded SQS entry. The HT-17 controlled fixtures exercised both successful and
intentional failure outcomes; they are state-machine evidence, not a claim of
scientific RS-LMTO convergence.

The installed-wheel legacy regression suite passed semantically for `rsutil
-S` and `rsutil -J`, covering ordinary bulk, multisite/multispecies, a
source-derived unequal-count lattice, seeded supercell/SQS, exchange
generation, and an explicit supported unique/SOC reduction. The reduction
assertion compares generated `njij` inventories and manifest context rather
than job-directory counts. Compatibility is workflow/semantic compatibility,
not byte-for-byte output identity.

### Provenance and limitations

No remote CI result was available from this checkout. The latest available
remote branch reference is `origin/HT` at
`54a279b6e00b5a51009626d94a16d1964149fc45`, which predates the closeout
commits; no remote pass is claimed. CI configuration for the source and
installed-wheel lanes is in `.github/workflows/build_pip.yaml`.

Pilot readiness is restored only within the support matrix above. HCP and
surface/impurity workflows remain restricted; full tensor exchange and
automatic ASE setup remain experimental; generic noncollinear symmetry,
conductivity, damping, AiiDA/UppASD integration, database backend, and ML
training remain unsupported. The earlier R09 and HT-18 readiness claims are
historical and superseded; they are preserved below for traceability.

## Historical HT-18 evidence (preserved)

> **REMEDIATION REQUIRED — THIS BRANCH IS NOT HT-READY FOR SCIENTIFIC
> PRODUCTION.**
>
> This document preserves the prior HT-18 closeout as historical evidence.
> The adversarial audit identified unresolved release blockers in source
> units, setup validation semantics, campaign state, execution evidence,
> manifest/reuse integrity, magnetic symmetry scope, source-format evidence,
> provenance, ASE setup policy, and clean-checkout reproducibility. See
> [`HT_REMEDIATION_BLUEPRINT.md`](HT_REMEDIATION_BLUEPRINT.md) for the frozen
> acceptance criteria and remediation ledger. The statuses and test results
> below must not be interpreted as a current production-readiness claim.

## Scope and provenance

This is the historical HT-18 closeout for the `HT` branch, retained for
baseline evidence only. The HT-18 work started from
the HT-17 tip `f40aa010710363c8a878f1d960c9b1b0799e5c15`; its immediate
parent is `be03295a953c4da8c94240962e83271dc293dddc`. The branch shares base
commit `e4a536d9660e3172ee5e6b4168eacf7b58331b4b` with `origin/main`.

The historical closeout was run on 2026-08-31 in the repository checkout with Python
3.11 and the existing compiled RS-utils library. The supplied untracked
`docs/dev/RSUTIL_HT_LUNA_BLUEPRINT_PACK/` directory was treated as campaign
input and was not modified. The protected Fortran sources under `source/`
were not changed.

## Historical HT-18 task ledger

| Task | Completed scope | Evidence |
| --- | --- | --- |
| HT-00 | Baseline and compatibility contract | `docs/HT_BASELINE_INVENTORY.md`, `docs/HT_COMPATIBILITY_CONTRACT.md` |
| HT-01 | Legacy semantic regression harness | `tests/test_legacy_cli.py` — 8 passed |
| HT-02 | Core typed contracts and safe paths | `tests/test_core_contracts.py` |
| HT-03 | Process-isolated Fortran setup runner | `tests/test_fortran_setup_runner.py` — 4 passed |
| HT-04 | Setup service behind `-S` | `tests/test_setup_service.py` |
| HT-05 | Transactional/idempotent setup | setup failure-injection tests |
| HT-06 | Validation, provenance, machine reports | `tests/test_setup_validation.py` |
| HT-07 | Canonical periodic `PairKey(i,j,R)` | `tests/test_ht07_pair_identity.py` |
| HT-08 | Exchange service behind `-J` | `tests/test_ht08_exchange_service.py` |
| HT-09 | Exact periodic and conservative SOC reduction | `tests/test_ht09_exchange_symmetry.py` |
| HT-10 | Transactional exchange jobs/manifests | `tests/test_ht10_exchange_jobs.py` |
| HT-11 | SCF/electronic result parser | `tests/test_ht11_results.py` |
| HT-12 | Exchange, DMI, and fail-closed tensor parsing | `tests/test_ht12_exchange_results.py` |
| HT-13 | Transport/damping contracts and unsupported path | `tests/test_ht13_transport_results.py` |
| HT-14 | Versioned TB model and ML export | `tests/test_ht14_tight_binding.py` |
| HT-15 | Optional ASE interoperability | `tests/test_ht15_ase.py` |
| HT-16 | Deterministic campaign generation | `tests/test_ht16_campaign.py` |
| HT-17 | Slurm rendering, status, and collection | `tests/test_ht17_campaign.py` |
| HT-18 | Readiness evidence and pilot closeout | This report and `examples/campaigns/ht18_pilot.json` |

## Architecture summary

The HT layer keeps the public `-S`/`-J` workflows and ordinary RS-LMTO file
names while moving reusable work behind path-explicit services. Setup inputs
are staged into an owned workspace, passed to a child Python worker, validated,
and atomically published. The legacy Fortran numerical kernel remains the
calculation authority; its fixed filenames and fatal-stop paths are contained
by the worker boundary.

Exchange identity is the typed periodic key `(i, j, R)`, with integer lattice
translation `R`; primitive sites are identified geometrically. Each site also
carries its source RS type from `izp`, while exchange-output endpoint integers
are RS type identities and are never interpreted as site numbers. Cartesian
displacement, distance, species, shell, and `Atom_*` directory names are
derived/provenance fields. Versioned result models preserve units, parser
provenance, raw source hashes, and validation
diagnostics. Campaign generation writes relocatable manifests and a stable
index; Slurm support only renders a reviewed array script and local status/
collection tools do not query or control a scheduler.

## Legacy compatibility evidence

The original HT-01 semantic corpus was rerun:

```text
python -m pytest -q tests/test_legacy_cli.py
8 passed in 30.15s
```

The corpus exercised real subprocess calls through `rsutil -S` and `rsutil
-J`, including simple magnetic setup, Fe-N multi-site setup, 1×1×2
supercell, seeded SQS, non-destructive/in-place behavior, Fe-Pt exchange
directories, Fermi-energy propagation, `Atom_*` layout, `--rij-out`, and the
non-crashing `--unique --soc` interface. Generated namelists were compared
semantically with `f90nml`; source inputs remained unchanged.

The broader legacy/worker subset also passed:

```text
python -m pytest -q tests/test_legacy_cli.py tests/test_fortran_setup_runner.py
12 passed in 33.70s
```

## Fortran numerical compatibility

The direct-versus-worker golden comparison used the checked-in bcc-Fe setup
fixture. The direct ctypes call and the child-worker call both succeeded. The
five generated files (`Fe.nml`, `input_new.nml`, `lattice_new.nml`,
`radii.out`, and `spheres.log`) matched byte-for-byte; the two generated
namelists also matched semantically. No Fortran source was modified.

The fatal-worker injection exited nonzero for an invalid element while the
parent test process continued and retained the isolated workspace. This is
covered by `test_fatal_worker_exit_is_returned_without_killing_parent`.

## Setup atomicity and failure injection

The targeted failure evidence passed 6 tests:

```text
python -m pytest -q \
  tests/test_fortran_setup_runner.py::test_fatal_worker_exit_is_returned_without_killing_parent \
  tests/test_setup_service.py::test_failed_worker_or_validation_does_not_publish_or_modify_source \
  tests/test_setup_service.py::test_staging_copy_failure_does_not_publish_or_modify_source \
  tests/test_ht10_exchange_jobs.py::test_unmarked_atom_directory_is_never_deleted_by_default \
  tests/test_ht12_exchange_results.py::test_mismatched_manifest_fails_instead_of_guessing
6 passed
```

The observed policies are:

* worker, validation, staging, and commit failures leave source namelists
  byte-identical and do not publish partial candidates;
* a stale/unmarked exchange directory is an explicit collision and is not
  deleted by default;
* a result row that cannot map uniquely to the supplied manifest/reference
  structure raises an interaction-mapping error rather than attaching to a
  guessed `PairKey`.

## Exchange pair and symmetry evidence

The exact symmetry/interop/campaign subset passed 32 tests (including pair
identity, exact symmetry, TB, ASE, and campaign coverage). Measured
symmetry values were:

| Fixture/check | Measurement |
| --- | --- |
| bcc Fe pair `PairKey(0,0,(1,0,0))` | 96 crystallographic operations; 12-member exact orbit; reverse translation present |
| B2 Fe-Co | 32 input pairs reduced to 2 representatives and expanded back to all 32 |
| Fe-Pt SOC reference | 48 non-SOC operations reduced to 8 operations preserving the +z axial-vector reference |
| MnSi scalar + DMI fixture | 584 pair-keyed records, `so` channel, validation passed |

Reverse semantics are explicitly `reverse(i,j,R) = (j,i,-R)`. Cartesian
displacements and species are derived from the key and structure; they are not
used as the identity.

## Parser coverage matrix

Statuses in this historical matrix mean `READY` was claimed by HT-18 based on
a checked-in stable fixture and tests;
`EXPERIMENTAL` is implemented but lacks complete upstream fixture evidence;
`DEFERRED` requires an upstream contract or future campaign work; and
`UNSUPPORTED` is deliberately fail-closed.

| Surface | Status | Evidence/limit |
| --- | --- | --- |
| SCF report/electronic records | READY | B2Fe and FePt reports parse with units, moments, provenance, and validation; convergence is only taken from explicit log markers |
| Scalar exchange (`jijso.out`/`jijfo.out`) | READY | Real MnSi fixture; canonical periodic mapping; audited source unit is mRy |
| DMI (`dijso.out`/`dijfo.out`) | READY | Real MnSi fixture; vector and SOC channel retained separately; audited source unit is mRy |
| Full exchange tensor | EXPERIMENTAL | `aij*.out` source scaling is audited as mRy; no checked-in fable_v4 tensor fixture; generic unverified tensors retain unknown units |
| Legacy `rs2txt.dat` reading | READY | Existing bands reader remains available; conversion to periodic TB is fail-closed when identity is ambiguous |
| Versioned fable_v4 TB | READY | `rs-tb-meta-v1`/`rs-tb-hoppings-v1` fixture parses and reconstructs complex `H(k)` |
| TB ML export | READY | Deterministic compressed NPZ/JSON; complex values and identity arrays round-trip |
| ASE structure adapter | READY | Optional boundary adapter preserves cell rows, PBC, species/order, metadata, and arrays when ASE is installed |
| Transport/damping text outputs | UNSUPPORTED | No evidenced stable RS-LMTO source/writer/header exists in this checkout; ambiguous fixture returns `UNSUPPORTED_AMBIGUOUS` without values |
| Generic noncollinear magnetic symmetry | DEFERRED | Only the documented collinear +z SOC reference is claimed |
| AiiDA/UppASD workflow integration | DEFERRED | Explicitly outside HT-18 scope |

## End-to-end generation/parse dry pilot

The reproducible pilot specification is
`examples/campaigns/ht18_pilot.json`. It contains four real setup entries:

| Entry | Coverage |
| --- | --- |
| `bcc-fe-reference` | simple elemental magnetic metal |
| `fen-multisite-reference` | multi-site, multi-species Fe-N structure |
| `fept-soc-reference` | SOC-relevant Fe-Pt setup/exchange metadata; FePt `control%nsp=2` is used for TB SOC evidence |
| `bcc-fe-seeded-supercell` | seeded lattice distortion with a 2×1×1 supercell |

Both independent generations completed with the real setup worker. They had
campaign ID
`b4ee6b0e78e5543ec08b8e9f63e88a1192b689e1cb1efba6e00a5011debc1e30`, spec
hash
`65a9de4bf7f242d5a938e0ba9a2e5e0cb023135640bd85de62012002ab1681a2`, and the
same entry IDs:

```text
bcc-fe-reference          7631e1a2efe058e6c866913d82a39d776c2306b30b060f3d116f7371190a8549
fen-multisite-reference   5b8e017b88d774fd41afb8620bcb69c51d2030aaa50f57c208f61614a87371bb
fept-soc-reference        406dfe8b1350b452d96faecdc0121710e108dc2eab7ef804fe6aabece33f7d14
bcc-fe-seeded-supercell   ad63b2d3a15cbb5804b40b4eae94f0c3fc11a94b3d509a3ffcc259c4ea231d6a
```

Across the two fresh destinations, all 23 non-manifest generated scientific
files had identical SHA-256 hashes. Run manifests intentionally contain
run-time timestamps and therefore are not expected to be byte-identical;
scientific identity and generated candidate content are deterministic.

The parser portion of the pilot exercised:

* B2Fe and FePt electronic reports: validation passed, two spin and two
  orbital moment records each;
* MnSi `jijso.out` + `dijso.out`: 584 merged canonical records, including
  `PairKey(0,0,(-2,0,0))` with `J=-0.00587 mRy` and
  `D=(0.000054,-0.000041,0.000138) mRy`;
* the versioned TB fixture: six hoppings, dimension four, eV units, complex
  `H(k)` reconstruction, `control%nsp=2`, SOC enabled, and lossless ML
  round-trip;
* the transport candidate fixture: explicit
  `UNSUPPORTED_AMBIGUOUS`, no guessed values, and a retained SHA-256 hash.

## Campaign and Slurm dry-run evidence

The historical four-entry index was rendered with a Python fake executable that
wrote only minimal `report.out` files. That evidence predates HT-R03: it did
not persist execution records and did not establish explicit SCF convergence.
Under the repaired state contract those entries are parser-valid at most, with
execution and scientific convergence unknown; they must not be described as
scientifically `COMPLETED_VALID`. The historical collection output therefore
does not constitute four successful scientific results.

The current HT-R03 regression matrix covers zero-exit converged,
non-converged, unknown, malformed, non-zero/no-output, stale-output, legacy
missing-record, and wrong-identity cases. The generated wrapper now persists
`rsutil_execution.json` atomically and collection exposes independent
execution, output/parser, and scientific axes. This validates the local
campaign boundary without claiming that Slurm was submitted or queried.

## Historical full verification

The full local suite passed:

```text
python -m pytest -q
97 passed, 10 warnings in 49.03s
```

The repository's existing CI workflows remain in `.github/workflows/`,
including `build.yaml` and `build_pip.yaml` (the latter runs `pytest -q`). No
remote CI run was available to this checkout; the local full-suite result is
historical closeout evidence only. It does not satisfy the current
clean-checkout/CI release gate.

## Known limitations and deferred work

* No fable_v4 source checkout is present here for transport/damping or full
  tensor writer evidence. Those parsers remain fail-closed or experimental.
* The legacy setup kernel can report `lld=21` for the default historical setup
  while the normalized spec carries `rec_steps=20`; this is retained as a
  warning for compatibility.
* The checked-in electronic reports contain explicit `Not converged!` log
  markers and the upstream report writer does not establish units for every
  possible energy/moment field. The parser does not infer convergence or
  units.
* TB overlap data is absent from the upstream versioned writer and remains
  absent in the model. `control%nsp` is recorded only when supplied by the
  input/metadata contract.
* Generic noncollinear magnetic symmetry, scheduler submission/retry/restart,
  database backends, ML training, AiiDA integration, and new UppASD workflow
  integration are not part of HT-18.

These limits are intentional. They should be resolved by adding upstream
format evidence and focused scientific tests before changing the support
matrix.
