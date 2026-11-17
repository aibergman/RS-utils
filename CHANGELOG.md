# Changelog

## 1.0.0 — 2026-11-17

First public release, distributed on PyPI as `rslmto-utils` under the GNU GPL v3. It follows the 0.11.0rc1 candidate (31 August 2026); the HT readiness evidence in `docs/HT_READINESS_REPORT.md` was gathered on that candidate. Changes since then:

### New workflows

- CIF to RS-LMTO conversion (replaces the external `cif2cell` step; see `docs/HT_CIF.md`):
  - CIF conversion now reduces centred cells to the standard primitive cell by default, matching `cif2cell` (e.g. 2 atoms for zinc-blende FeN instead of the 8-atom cubic cell). Use `--cell conventional` to keep the CIF cell; this replaces the `--primitive` flag of the unreleased converter.
  - `rsutil setup cif` / `rsutil.cif`: convert ordered CIF files to RS-LMTO `lattice.nml` and `input.nml` (pymatgen-based, one RS type per atom), validated against `cif2cell` output. See `docs/HT_CIF.md`.
  - `rsutil setup bulk --cif FILE` converts a CIF into the working directory and then runs the normal setup; `ct` is chosen automatically from a target neighbour count (`--ct-neighbors`, default 15) unless `--ct` is given.
- Surface preparation from a converged bulk calculation (`rsutil setup surface`): Miller-index surfaces, layers and empty-sphere types, `--orientation global-z`, in-plane bond alignment, and expansion of an existing surface.
- Embedded impurity setup in bulk and surface hosts (`rsutil setup impurity`).
- Campaign exchange execution stage (`rsutil campaign exchange prepare|render-slurm|status|collect`), campaign scaffolding from collections of RS-LMTO source bundles, heterogeneous per-entry Slurm resources, and scalar summaries.
- DOS plotting from calculation directories (`rsutil viz dos`) with interactive controls.
- Single-run workflows are grouped under `rsutil setup`; `-S` remains as a shortcut.

### Setup, exchange and export

- `--preset KS` setup option and a default `strux_backend`; fixes to first-pass empty-sphere generation and to `ct`/`r2` allocation in the Fortran setup kernel.
- Exchange pair slicing and magnetic defaults; UppASD `momfile` moments are read from `report.out`.
- CIF provenance is kept separate from campaign entry identity.

### Packaging, licensing and robustness

- Packaging for public release: distribution renamed `rslmto-utils` (import name and command remain `rsutil`), metadata moved to `pyproject.toml`, GPL-3.0 `LICENSE` and a `NOTICE` recording the provenance of derived Fortran sources, minimum Python raised to 3.11, clear error when CMake or a Fortran compiler is missing for source builds.
- Failures in band-analysis numerics (Chern numbers, effective masses, tetrahedra) now warn and yield NaN instead of being silently replaced by plausible values; removed a dead call to a non-existent Berry-curvature method.
- Importing `rsutil` no longer warns when the Fortran library is missing; the 21 kernel tests are skipped, not failed, without it.

## 0.11.0rc1 — 2026-08-31

Release candidate for controlled high-throughput (HT) RS-LMTO pilot use.

### High-throughput foundation

- Transactional, isolated `-S` setup with a protected Fortran worker.
- Typed structure, setup, and result contracts with deterministic manifests and provenance.
- Lightweight campaign generation, filesystem status, collection, and Slurm array rendering.
- Convergence-aware run-state handling that keeps execution, parsing, and scientific convergence distinct.

### Exchange

- Canonical periodic `PairKey(i,j,R)` identity and RS site-to-RS-type preservation through `izp`.
- Scalar exchange `J`, DMI, anisotropic `Aij`, and derived full interaction tensor contracts.
- Explicit source units, conservative magnetic symmetry reduction, and hardened exchange-job manifests/reuse.

### Data and interoperability

- Versioned tight-binding/ML data formats and ASE `Atoms <-> RSStructure` interoperability.
- Machine-readable JSON/NPZ results and explicit unsupported or restricted domains for formats without a safe contract.

### Compatibility and qualification

- Traditional `rsutil -S` and `rsutil -J` workflows remain supported.
- `0.11.0rc1` is intended for controlled HT pilot use. Some workflows remain `RESTRICTED`, `EXPERIMENTAL`, or `UNSUPPORTED` according to the readiness matrix; this release does not promise universal RS-LMTO workflow coverage.
