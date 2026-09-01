# RS-LMTO transport and damping format audit

## Audit boundary

This document records the HT-R06 audit of the public `fable_v4` branch of
`rslmtoasa/rslmtoasa`. The branch was fetched on 2026-08-31 and pinned for
this audit at commit
`faf71320c14afd680ce6973fb48fe4b6f8a25529`. The source links below retain
the branch name so a later audit can compare the pinned revision with the
current public branch.

The result is deliberately conservative. `fable_v4` does contain transport
and damping writers; the previous claim that they were absent was incorrect.
However, the emitted text does not carry enough machine-readable units and
run context for RS-utils to parse it safely as an HT result. No Python
transport physics or unit conversion is added here.

Audited source references:

* [`conductivity.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/fable_v4/source/conductivity.f90#L478-L766)
  — conductivity reconstruction and all output writers;
* [`recursion_transport.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/fable_v4/source/recursion_transport.f90#L368-L890)
  — KPM moments and charge/spin/orbital/accumulation/torque operator modes;
* [`control.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/fable_v4/source/control.f90#L205-L223)
  and [`control` namelist declarations](https://github.com/rslmtoasa/rslmtoasa/blob/fable_v4/source/include_codes/namelists/control.f90#L1-L23)
  — transport controls and defaults;
* [`calculation.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/fable_v4/source/calculation.f90#L754-L852)
  and [`calculation` namelist declarations](https://github.com/rslmtoasa/rslmtoasa/blob/fable_v4/source/include_codes/namelists/calculation.f90#L1-L9)
  — damping dispatch and gates;
* [`exchange_dynamics.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/fable_v4/source/exchange_dynamics.f90#L1-L282)
  — Gilbert and inertia kernels/writers;
* [`spin_dynamics.f90`](https://github.com/rslmtoasa/rslmtoasa/blob/fable_v4/source/spin_dynamics.f90#L65-L90)
  — phenomenological LLG `alpha` input.

## Source-grounded inventory

| Observable / output | Status | Source, controls, and emitted schema | Why this is or is not HT-parseable |
| --- | --- | --- | --- |
| Charge KPM/Kubo response: `cond_total.out` | `UNSUPPORTED_AMBIGUOUS` | `source/conductivity.f90::calculate_conductivity_tensor` declares the file and writes three fields: `E-Ef`, integrated real response, integrated imaginary response. It is reached with `&calculation post_processing='conductivity'`; `&control cond_type='charge'`, `cond_ll`, `cond_calctype`, and `random_vec_num` control the route. | The writer has no header, unit, method/options record, broadening/regularization record, volume/unit-cell convention, or tensor component labels. The active prefactor is `16/(pi*de**2)`; the volume/electronic prefactor alternatives are comments, not an emitted contract. A three-column numeric file therefore cannot safely be labeled electrical conductivity or assigned SI/code units. |
| Spin KPM/Kubo response: `cond_total.out`, `<element>_cond.out` | `UNSUPPORTED_AMBIGUOUS` | `source/recursion_transport.f90::compute_moments_stochastic` selects `cond_type='spin'` and substitutes `js_a/jso_a` for the charge velocity operators. `cond_calctype='per_type'` emits `<element>_cond.out`; `random_vec` uses `random_vec_num`. | The output shape remains the same three columns as charge and does not identify spin-current convention, spin axis, normalization, or units. The filename and `cond_type` are not present in the file, so a parser cannot recover the channel from output bytes alone. |
| Orbital KPM/Kubo response: `cond_total.out`, `<element>_cond.out` | `UNSUPPORTED_AMBIGUOUS` | The same recursion routine selects `cond_type='orbital'` and substitutes `jl_a/jlo_a`. `jl_alpha` selects `x`, `y`, or `z` for accumulation/operator paths; the normal conductivity output still has no channel metadata. | The writer does not preserve the orbital-current convention, angular-momentum normalization, component order, or units. It is unsafe to treat this as a conductivity result merely because it is written by `conductivity.f90`. |
| Orbital-resolved real/imag response: `cond_total_orb_real.out`, `cond_total_orb_im.out`, and per-type `<element>_cond_orb_real.out` / `<element>_cond_orb_im.out` | `UNSUPPORTED_AMBIGUOUS` | `calculate_conductivity_tensor` writes `E-Ef` followed by `integrated_l_real(:,i)` or `integrated_l_im(:,i)` using `(19es16.6)`. The code's `l2=1..nb` loop supplies the component sequence; per-type records use `integrated_type_l_*`. | There is no header or `lmax`/basis map in the file. The fixed 19-field format is not a self-describing component contract, and the output still lacks units and operator/channel metadata. No component transpose or orbital parser is justified. |
| Spin/orbital accumulation and torque response modes | `UNSUPPORTED_AMBIGUOUS` | `compute_moments_stochastic` also accepts `spin_accumulation`, `orbital_accumulation`, `spin_torque`, and `orbital_torque`, using `js_alpha`, `jl_alpha`, and the corresponding operator builders. They reuse the conductivity output writer and filenames. | These are response/operator modes, not established electrical conductivity formats. Their output has the same unlabeled numeric schema, with no response observable name, axis, units, or normalization. They must not be promoted into `ConductivityResult`. |
| Frequency-dependent / optical conductivity with an explicit `omega` axis | `NOT_FOUND` | The current source has `energy_min`, `energy_max`, `channels_ldos`, and emits `this%en%ene(nv)-this%en%fermi`; no frequency grid or `omega` field is consumed by the writer. The narrative reference page describes optical conductivity, but the source writer does not establish that output contract. | Do not reinterpret the energy-resolved cumulative real/imaginary files as an optical-frequency spectrum. A future frequency writer needs an explicit frequency axis and units. |
| Pair Gilbert/Kamberský damping: `alldampings.out` | `UNSUPPORTED_AMBIGUOUS` | `source/exchange_dynamics.f90::calculate_gilbert_damping` is called from `source/calculation.f90::post_processing_exchange` when `&calculation do_damping=.true.`. It is gated by `nsp=2` and source comments require SOC; k-space routes additionally require second-order `H(k)`. The header is `i j`, nine values in row-major nested `k=1..3`, `l=1..3` order, `0.5*(xx+yy)`, distance, and three Cartesian `rij` values. The pair indices are the source lattice's one-based indices. | The writer provides useful column labels but no units, g-factor/normalization record, SOC parameters, Green-function broadening/route, lattice-length unit for distance/`rij`, or input/run identity. The source calls the values Gilbert damping, but a filename/header alone cannot establish the complete quantitative contract. |
| Energy-resolved total damping: `damping-energy.out` | `UNSUPPORTED_AMBIGUOUS` | The same subroutine writes `E-Ef` plus nine real tensor components after summing pair kernels. The values are sampled on the electronic energy mesh; `ief` selects the nearest Fermi point for the pair table. | No units or tensor normalization are emitted. The active source uses a single `factor` when writing the accumulated energy table after the pair loop, so the total table's normalization is not independently recoverable from the file for heterogeneous pair moments. |
| Experimental magnetic inertia kernel: `inertia-energy.out`, `allinertias.out`, `example-real.out`, `example-imag.out` | `UNSUPPORTED_AMBIGUOUS` | `source/exchange_dynamics.f90::calculate_moment_of_inertia` writes these files under `do_inertia=.true.`. Its own source comment calls this an experimental raw torque-correlation kernel with no independently established production relation or physical normalization. | These are not Gilbert damping. They have raw real/imaginary derivative-response columns and no established units/normalization; RS-utils must not expose them as `DampingResult`. |
| Atomistic LLG damping input / trajectory: `spin_dynamics.f90`, `output.lammpstrj` | `NOT_FOUND` as an ab-initio damping output | `spin_dynamics.f90` reads phenomenological `alpha` from its spin-dynamics input and writes a trajectory. This is a time-integration input/trajectory path, not a writer for an ab-initio Gilbert tensor. | A relaxation parameter or trajectory must not be reinterpreted as calculated Gilbert damping. |

## Evidence details

### Conductivity and KPM route

The `conductivity` type reads only `ll_cond` from the separate
`&conductivity` namelist, while the operational controls are in `&control`:
`cond_ll`, `cond_type`, `cond_calctype`, and `random_vec_num`. The source
defaults `cond_ll=200`, `cond_type='charge'`, and `cond_calctype='per_type'`.
`recursion_transport.f90` constructs two-index stochastic Chebyshev moments
and documents the route as charge, spin, orbital, accumulation, and torque
operator insertions.

`calculate_conductivity_tensor` uses an energy mesh of
`channels_ldos + 10`, rescales it using `energy_min`/`energy_max`, integrates
the real and imaginary arrays up to each Fermi-weighted energy point, and
writes plain formatted records. The total and per-type records are
`E-Ef, real, imaginary`; orbital records are `E-Ef` plus the `l2` sequence.
The source does not write a header or any of the required HT interpretation
metadata. In particular, `volume` is computed but the active `factor` omits
it; the alternative expressions involving physical constants are commented
out. This is sufficient to reproduce the source's numeric output, not to
assign a portable conductivity unit.

The source also has no independent `omega`, temperature, Fermi-Dirac, or
finite-frequency output axis in this writer. The Lorentz kernel parameter
`6.0d0` and Chebyshev order are implementation details, but are not emitted
with each result. The current RS-utils `ConductivityResult` contract requires
all such interpretation metadata, so no legacy text parser is registered.

### Gilbert damping and damping-adjacent paths

`post_processing_exchange` routes `do_damping` to
`calculate_gilbert_damping` after the selected recursion/Lehmann/Dyson Green
function has filled the canonical intersite arrays. The damping routine calls
`torque_operator_collinear`, forms anti-Hermitian Green-function blocks, takes
real and imaginary traces for the nine `(k,l)` tensor components, and applies
`(-0.25)*(2/(pi*spin_i))` at the nearest Fermi-energy grid point for
`alldampings.out`. These are concrete source semantics, not an absence claim.

They still do not make a safe result format: the text files do not identify
the Green-function route or broadening, do not preserve SOC and torque
operator provenance, and do not declare the units of tensor, distance, or
displacement fields. `calculate_moment_of_inertia` explicitly warns that its
kernel is experimental and raw. The phenomenological LLG `alpha` input and
TDDFT/response imaginary parts are separate concepts and are not promoted to
Gilbert damping here.

## RS-utils status and tests

`rsutil.results.transport` retains typed `ConductivityResult` and
`DampingResult` contracts for a future evidenced format. Its registry now
records the fable_v4 source evidence while remaining
`UNSUPPORTED_AMBIGUOUS`. The parser recognizes the actual fable_v4 filenames
so a directory containing `cond_total.out`, `<element>_cond.out`,
`damping-energy.out`, or `alldampings.out` returns an explicit unsupported
result with hashes and no guessed values. Existing generic candidate names
retain the same fail-closed behavior.

No `SUPPORTED_STABLE` or `SUPPORTED_EXPERIMENTAL` parser is claimed. No real
or source-faithful fixture is promoted because a fixture cannot repair the
missing units and metadata in the writer contract.

## Proposed upstream machine-readable contract

For any future HT parser, fable_v4 should emit a versioned sidecar or record
containing at least:

* observable and format version;
* run identity, input/options hash, Green-function route, and source revision;
* method and approximation, including Chebyshev order, kernel/broadening,
  temperature, and Fermi convention;
* every energy or frequency axis with its unit and whether it is `E-Ef`, `E`,
  or `omega`;
* conductivity channel (`charge`, `spin`, `orbital`, accumulation, or
  torque), operator-axis choices, component order, and units;
* per-type/site index convention and an explicit basis/orbital map;
* damping observable name, site/pair identity convention, tensor order,
  normalization, SOC/torque provenance, and units;
* an explicit distinction between Gilbert damping, magnetic inertia,
  relaxation rates, and generic imaginary response.

Until that contract exists, RS-utils remains fail-closed.
