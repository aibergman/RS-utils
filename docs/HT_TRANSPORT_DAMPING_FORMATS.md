# HT-13 transport and damping formats

## Evidence boundary

The checked-out RS-utils repository contains no checked-out `fable_v4` source
tree. Its local `source/` tree is the setup/sphere/potential-parameter kernel;
the complete inventory contains no `conductivity.f90`, Kubo/KPM module,
transport writer, Gilbert-damping module, or damping writer. The checked-in
examples likewise contain no transport or damping output fixture.

Consequently, HT-13 declares no production parser `SUPPORTED_STABLE`. This is
intentional: a numeric table without an evidenced method, channel, axes, units,
and component/index ordering must not become a confident high-throughput
result.

## Classification

| Observable/output surface | Status | Evidence and reason |
| --- | --- | --- |
| `conductivity.f90` / Kubo output | `UNSUPPORTED_AMBIGUOUS` | No such source or output header is present in this checkout. The required upstream method, broadening, channel, tensor order, and units are unavailable. |
| KPM transport output | `UNSUPPORTED_AMBIGUOUS` | No KPM source, namelist option, report header, or fixture is present. A filename containing `kpm` is not sufficient evidence of a schema. |
| Gilbert/magnetic damping output | `UNSUPPORTED_AMBIGUOUS` | No damping implementation or output path is present. Site/pair indexing, tensor order, normalization, and the distinction between a Gilbert tensor and another response are therefore unknown. |
| `report.out` electronic records | `UNSUPPORTED_AMBIGUOUS` for HT transport/damping | The existing stable parser is intentionally limited to electronic/SCF records such as Fermi energy and moments. It contains no transport/damping observable and is not reinterpreted here. |
| `rsutil.bands` Berry/spin conductivity helpers | `UNSUPPORTED_AMBIGUOUS` as an RS-LMTO result format | These are Python post-processing helpers, not a documented `fable_v4` output contract. Their approximations and units cannot establish a parser for RS-LMTO transport files. |

No `SUPPORTED_EXPERIMENTAL` parser is registered either. There is not enough
format evidence to implement even a version-gated text parser safely.

## Implemented contract

`rsutil.results.transport` provides a small, non-plugin registry and typed
future-facing models:

* `ConductivityResult` requires explicit result axes, a unit-bearing value,
  component order, method, channel, approximation, options, and source
  provenance. Its value shape is checked as `axes + component_order`.
* `DampingResult` requires explicit axes, `site` or `pair` indexing, index
  labels, component order, method, channel, approximation, options, and
  provenance. Its value shape is checked as `axes + index_labels +
  component_order`.
* `parse_transport_result()` and its conductivity/damping aliases return an
  `UnsupportedTransportResult` with status `UNSUPPORTED_AMBIGUOUS` for a
  candidate file. They never return partial guessed values. Missing explicit
  files still raise the normal `MissingResultFileError`.
* Provenance includes SHA-256 hashes for the candidate output and neighboring
  `input.nml`, `lattice.nml`, and manifest files when present.

The registry intentionally has zero stable parsers at this point, so there is
no claimed stable fixture. The tests include an ambiguous transport fixture to
verify the explicit unsupported path and provenance hashes.

## Minimal future `fable_v4` change

The upstream implementation should emit one versioned machine-readable result
record (JSON, CSV plus a sidecar, or equivalent) containing at least:

* observable and format version;
* method, approximation, broadening/regularization, and spin/SOC channel;
* every axis and its unit (energy/Fermi grid where applicable);
* conductivity component order and units;
* damping site/pair index convention, tensor component order, and units;
* run identity and input/options hash.

That change belongs in the RS-LMTO repository, not in this RS-utils task.
