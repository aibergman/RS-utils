# Changelog

## Unreleased

- `rsutil setup cif` / `rsutil.cif`: convert ordered CIF files to RS-LMTO `lattice.nml` and `input.nml` (pymatgen-based, one RS type per atom), validated against `cif2cell` output. See `docs/HT_CIF.md`.
- `rsutil setup bulk --cif FILE` converts a CIF into the working directory and then runs the normal setup; `ct` is chosen automatically from a target neighbour count (`--ct-neighbors`, default 15) unless `--ct` is given.

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
