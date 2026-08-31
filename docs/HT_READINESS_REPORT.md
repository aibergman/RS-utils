# RS-utils HT readiness report — historical HT-18 evidence

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
translation `R`. Cartesian displacement, distance, species, shell, and
`Atom_*` directory names are derived/provenance fields. Versioned result
models preserve units, parser provenance, raw source hashes, and validation
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
