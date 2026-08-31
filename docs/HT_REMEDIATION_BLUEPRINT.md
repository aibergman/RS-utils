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
- [ ] **HT-R02** — Restore complete RS-LMTO setup/lattice validation semantics
- [ ] **HT-R03** — Implement scientific campaign state and execution evidence
- [ ] **HT-R04** — Harden manifests, pair identity, and reuse integrity
- [ ] **HT-R05** — Bound and harden exchange symmetry/magnetic-reference semantics
- [ ] **HT-R06** — Complete transport/damping source audit against `fable_v4`
- [ ] **HT-R07** — Provenance, ASE, path, indexing, and dead-code cleanup
- [ ] **HT-R08** — Make packaging, test environment, and CI reproducible
- [ ] **HT-R09** — True adversarial HT closeout and pilot evidence

No task may be marked complete before its acceptance tests pass. R02-R09 are
intentionally unchecked at the completion of HT-R01.

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
