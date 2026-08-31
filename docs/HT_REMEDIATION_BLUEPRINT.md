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
- [ ] **HT-R05** — Bound and harden exchange symmetry/magnetic-reference semantics
- [ ] **HT-R06** — Complete transport/damping source audit against `fable_v4`
- [ ] **HT-R07** — Provenance, ASE, path, indexing, and dead-code cleanup
- [ ] **HT-R08** — Make packaging, test environment, and CI reproducible
- [ ] **HT-R09** — True adversarial HT closeout and pilot evidence

No task may be marked complete before its acceptance tests pass. R04-R09 remain
intentionally unchecked until their respective repairs are completed.

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
