# Questaal empty-sphere finder

This directory is the vendored source for the optional Questaal empty-sphere
backend. `findes.f` is copied verbatim from Questaal's tracked
`src/subs/findes.f` at repository revision `b5c8b3cad3d2afe8fc71915c0c4231fabd7e9a20`
(blob `b6934c18b0830e32bdb2b22b5869f8bb1f33728c`).

The standalone radius-resizing adapter is a source-level adaptation of
Questaal's `src/subs/sclwsr.f` at the same revision (blob
`3216144d314fa37ec9fefd592d60b56e425cfcaa`). It retains the `sclws2`
species-locking and constrained-volume iteration while replacing the runtime
`pairc` neighbor-table builder with an in-library periodic-image enumerator.
The `nrmliz` shim follows `src/subs/mkplat.f` at the same revision (blob
`00bf58d45fcb9bfa184ffb37ce3e2621763897f9`).
The compatibility bridge retains FINDES' requested packing fraction across
the resize call so its outer loop can add spheres when that fraction is not
reachable with the current sites. The original `sclwsr` reports its achieved
fraction through the same argument, which would make FINDES' next target
comparison tautological. This is an adapter fix at the call boundary; the
pinned `findes.f` source remains byte-for-byte intact. The RS-utils adapter
sets the overlap limits to 0.14 of center distance and 0.26 of sphere radius;
zero limits caused FINDES to stop at tangent spheres and leave the FeN
XBand-style sphere-volume target underfilled.

The source file documents that the finder was adapted from the Stuttgart
package and records Questaal updates from 2011 through 2017. The Questaal
repository's `COPYING` applies AGPL-3.0-or-later terms. Preserve that notice;
RS-utils now declares AGPL-3.0-or-later in its root `LICENSE`.

The free-form copies in local Questaal working directories are not the source
of record; this snapshot is taken from the revisioned fixed-form source. The
adapter work should be recorded separately from this upstream snapshot.

## Extraction boundary

The public RS-utils adapter calls the finder in memory and returns the
expanded basis and sphere radii. It does not invoke `sufindes` or the
Questaal file-writing path. The finder source itself contains the geometric
helpers (`rsmesh`, `mdeses`, `fillat`, `bigges`, `symes`, and related routines),
but the top-level workflow also reaches RS-utils-external services:

- `sclwsr` performs radius resizing and overlap control and depends on
  Questaal's pair/neighbor infrastructure;
- `structures` supplies `str_ctrl` and `str_lat` used by the legacy entry;
- shortening, basis ordering, diagnostics, and symmetry-variable routines are
  supplied by the Questaal runtime.

The standalone adapter therefore uses explicit host-side structure types and
preserves the `sclws2` iterative radius/overlap policy while replacing the
`pairc` neighbor-table service. It receives the identity symmetry operation,
matching RS-utils' current xband symmetry handling; general symmetry
operations are not yet exposed by the RS-utils data model.
