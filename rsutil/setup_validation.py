"""Semantic setup validation and versioned provenance reports.

The validator deliberately works on the namelists rather than on console
output.  It is conservative about legacy inputs: fields that the Fortran
kernel supplies by default remain optional, while contradictions in fields
that are present are reported as coded diagnostics.
"""

from __future__ import annotations

import json
import math
import platform
import subprocess
import sys
import sysconfig
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import f90nml

from .core.contracts import (
    Diagnostic,
    DiagnosticSeverity,
    ErrorCode,
    OperationStatus,
    RunManifest,
    SetupReport,
    SetupResult,
    SetupSpec,
    ValidationResult,
)
from .core.hashing import hash_file, scientific_hash


_LATTICE_REQUIRED = (
    "nbulk_bulk",
    "ntot",
    "nbas",
    "nrec",
    "a",
    "crd",
    "izp",
    "no",
    "iu",
    "ib",
    "irec",
)
_KNOWN_ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
    "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
    "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
_ATOMIC_NUMBERS = {symbol: number for number, symbol in enumerate(_KNOWN_ELEMENTS, start=1)}


def _diag(
    severity: DiagnosticSeverity,
    code: ErrorCode | str,
    message: str,
    *,
    path: str | None = None,
    **details: Any,
) -> Diagnostic:
    return Diagnostic(severity, code, message, path=path, details=details)


def _error(code: ErrorCode, message: str, *, path: str | None = None, **details: Any) -> Diagnostic:
    return _diag(DiagnosticSeverity.ERROR, code, message, path=path, **details)


def _warning(code: ErrorCode, message: str, *, path: str | None = None, **details: Any) -> Diagnostic:
    return _diag(DiagnosticSeverity.WARNING, code, message, path=path, **details)


def _info(code: ErrorCode, message: str, *, path: str | None = None, **details: Any) -> Diagnostic:
    return _diag(DiagnosticSeverity.INFO, code, message, path=path, **details)


def _has_errors(diagnostics: Sequence[Diagnostic]) -> bool:
    return any(item.severity is DiagnosticSeverity.ERROR for item in diagnostics)


def _stochastic_provenance(spec: SetupSpec) -> dict[str, Any]:
    """Describe only stochastic operations active in ``spec``.

    A seed is not required for deterministic setup work.  The legacy seed
    applies to lattice/basis distortions, while SQS has its own seed.
    """

    operations: list[str] = []
    seeds: dict[str, int | None] = {}
    if spec.lattice_dist > 0:
        operations.append("lattice_distortion")
        seeds["lattice_distortion"] = spec.seed
    if spec.basis_dist > 0:
        operations.append("basis_distortion")
        seeds["basis_distortion"] = spec.seed
    if spec.sqs is not None:
        operations.append("sqs")
        seeds["sqs"] = spec.sqs_seed

    unseeded = [name for name in operations if seeds[name] is None]
    reproducible = not unseeded
    if not operations:
        reason = "no stochastic operations requested"
    elif reproducible:
        reason = "all active stochastic operations are seeded"
    else:
        reason = "unseeded stochastic operations: " + ", ".join(unseeded)
    return {
        "active": bool(operations),
        "operations": operations,
        "seeds": seeds,
        "seed": spec.seed,
        "sqs_seed": spec.sqs_seed,
        "reproducible": reproducible,
        "reason": reason,
    }


def _read_namelist(path: Path, diagnostics: list[Diagnostic]) -> Mapping[str, Any] | None:
    if not path.is_file():
        diagnostics.append(_error(ErrorCode.MISSING_INPUT, "required namelist is missing", path=str(path)))
        return None
    try:
        return f90nml.read(path)
    except Exception as exc:
        diagnostics.append(
            _error(
                ErrorCode.NAMELIST_PARSE_ERROR,
                f"could not parse namelist: {exc}",
                path=str(path),
            )
        )
        return None


def _group(
    data: Mapping[str, Any], name: str, diagnostics: list[Diagnostic], path: Path
) -> Mapping[str, Any] | None:
    for key, value in data.items():
        if str(key).lower() == name.lower() and isinstance(value, Mapping):
            return value
    diagnostics.append(
        _error(
            ErrorCode.NAMELIST_GROUP_MISSING,
            f"required namelist group '&{name}' is missing",
            path=str(path),
            group=name,
        )
    )
    return None


def _value(group: Mapping[str, Any], name: str) -> Any:
    for key, value in group.items():
        if str(key).lower() == name.lower():
            return value
    return None


def _array(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _as_int(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        converted = int(value)
    except (TypeError, ValueError):
        return None
    if isinstance(value, str):
        try:
            return converted if float(value.strip()) == converted else None
        except ValueError:
            return None
    return converted if converted == value else None


def _finite_number(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _determinant(vectors: Sequence[Sequence[Any]]) -> float | None:
    if len(vectors) != 3 or any(len(vector) != 3 for vector in vectors):
        return None
    try:
        a, b, c = (tuple(float(item) for item in vector) for vector in vectors)
    except (TypeError, ValueError):
        return None
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _label_element(label: str) -> tuple[str | None, int | None]:
    text = str(label).strip()
    if not text:
        return None, None
    if len(text) >= 2 and text[:2].capitalize() in _ATOMIC_NUMBERS:
        symbol = text[:2].capitalize()
    elif text[:1].capitalize() in _ATOMIC_NUMBERS:
        symbol = text[:1].capitalize()
    else:
        return None, None
    return symbol, _ATOMIC_NUMBERS[symbol]


def _labels(value: Any) -> list[str]:
    return [str(item).strip() for item in _array(value)]


def _validate_lattice_group(
    group: Mapping[str, Any],
    path: Path,
    diagnostics: list[Diagnostic],
    *,
    calculation_type: str | None = None,
    nlay: int | None = None,
    input_ntype: int | None = None,
) -> dict[str, Any]:
    """Validate lattice fields using their independent RS-LMTO meanings.

    The four count fields are deliberately not interchangeable.  ``ntot``
    sizes primitive-cell tables, ``nrec`` sizes recursion representatives,
    ``nbas`` is a charge-distribution/shell count, and ``nbulk_bulk`` is a
    mode input used by surface/interface construction.
    """
    normalized: dict[str, Any] = {}
    for name in _LATTICE_REQUIRED:
        value = _value(group, name)
        if value is None:
            diagnostics.append(
                _error(
                    ErrorCode.NAMELIST_FIELD_MISSING,
                    f"required lattice field '{name}' is missing",
                    path=str(path),
                    field=name,
                )
            )
        else:
            normalized[name] = value

    counts: dict[str, int] = {}
    for name in ("nbulk_bulk", "ntot", "nbas", "nrec"):
        value = _value(group, name)
        converted = _as_int(value)
        if converted is None or converted <= 0:
            diagnostics.append(
                _error(
                    ErrorCode.COUNT_MISMATCH,
                    f"lattice field '{name}' must be a positive integer",
                    path=str(path),
                    field=name,
                    value=value,
                )
            )
        else:
            counts[name] = converted
            normalized[name] = converted
    optional_nbulk = _value(group, "nbulk")
    nbulk = _as_int(optional_nbulk) if optional_nbulk is not None else None
    if optional_nbulk is not None:
        if nbulk is None or nbulk < 0:
            diagnostics.append(
                _error(
                    ErrorCode.COUNT_MISMATCH,
                    "lattice field 'nbulk' must be a non-negative integer",
                    path=str(path),
                    field="nbulk",
                    value=optional_nbulk,
                )
            )
        else:
            normalized["nbulk"] = nbulk

    # In the source, these arrays have distinct owners.  ``crd`` is the
    # primitive-cell basis, ``iu``/``ib`` identify the primitive-cell sites
    # used for bulk/reference calculations, and ``irec`` identifies the
    # recursion representatives.  Cluster indices cannot be upper-bounded
    # here because the cluster is built only after preflight validation.
    primitive_count = counts.get("ntot")
    recursion_count = counts.get("nrec")

    vectors = _array(_value(group, "a"))
    if len(vectors) != 3 or any(not isinstance(vector, Sequence) or len(vector) != 3 for vector in vectors):
        diagnostics.append(
            _error(
                ErrorCode.LATTICE_NONFINITE,
                "lattice field 'a' must contain three finite 3-vectors",
                path=str(path),
            )
        )
    elif any(not _finite_number(item) for vector in vectors for item in vector):
        diagnostics.append(
            _error(ErrorCode.LATTICE_NONFINITE, "lattice vectors must contain only finite numbers", path=str(path))
        )
    else:
        determinant = _determinant(vectors)
        if determinant is None or not math.isfinite(determinant):
            diagnostics.append(_error(ErrorCode.LATTICE_NONFINITE, "lattice determinant is not finite", path=str(path)))
        elif abs(determinant) <= 1.0e-10:
            diagnostics.append(
                _error(
                    ErrorCode.LATTICE_SINGULAR,
                    "lattice vectors are singular or near-singular",
                    path=str(path),
                    determinant=determinant,
                )
            )
        normalized["a"] = [list(vector) for vector in vectors]

    coordinates = _array(_value(group, "crd"))
    if primitive_count is not None and len(coordinates) != primitive_count:
        diagnostics.append(
            _error(
                ErrorCode.COUNT_MISMATCH,
                f"coordinate count {len(coordinates)} does not match ntot={primitive_count}",
                path=str(path),
                field="crd",
                expected=primitive_count,
                actual=len(coordinates),
            )
        )
    if any(not isinstance(vector, Sequence) or len(vector) != 3 for vector in coordinates):
        diagnostics.append(_error(ErrorCode.LATTICE_NONFINITE, "coordinates must contain 3-vectors", path=str(path)))
    elif any(not _finite_number(item) for vector in coordinates for item in vector):
        diagnostics.append(_error(ErrorCode.LATTICE_NONFINITE, "coordinates must be finite", path=str(path)))
    else:
        normalized["crd"] = [list(vector) for vector in coordinates]

    for name, expected in (
        ("izp", primitive_count),
        ("no", primitive_count),
        ("iu", primitive_count),
        ("ib", primitive_count),
        ("irec", recursion_count),
    ):
        values = _array(_value(group, name))
        if expected is not None and len(values) != expected:
            diagnostics.append(
                _error(
                    ErrorCode.COUNT_MISMATCH,
                    f"array '{name}' has {len(values)} entries; expected {expected}",
                    path=str(path),
                    field=name,
                    expected=expected,
                    actual=len(values),
                )
            )
    for name in ("izp", "no", "iu", "ib", "irec"):
        values = _array(_value(group, name))
        if values and any(_as_int(value) is None or int(value) < 1 for value in values):
            diagnostics.append(
                _error(
                    ErrorCode.COUNT_MISMATCH,
                    f"array '{name}' must contain positive integer indices",
                    path=str(path),
                    field=name,
                )
            )

    if _value(group, "ntype") is not None:
        ntype = _as_int(_value(group, "ntype"))
        if ntype is None or ntype <= 0:
            diagnostics.append(_error(ErrorCode.COUNT_MISMATCH, "lattice ntype must be positive", path=str(path)))
        else:
            normalized["ntype"] = ntype

    # ``ntype = nbulk + nrec`` is a source relationship when nbulk is part of
    # the materialized lattice state.  File-mode inputs often omit nbulk, so
    # do not manufacture it merely to make this check possible.
    if nbulk is not None and input_ntype is not None and recursion_count is not None:
        expected_ntype = nbulk + recursion_count
        if input_ntype != expected_ntype:
            diagnostics.append(
                _error(
                    ErrorCode.COUNT_MISMATCH,
                    f"input ntype={input_ntype} does not match source nbulk+nrec={expected_ntype}",
                    path=str(path),
                    field="ntype",
                    expected=expected_ntype,
                    actual=input_ntype,
                    mode=calculation_type,
                )
            )

    # Surface/interface construction derives the active recursion count from
    # nlay.  Apply these checks only when the mode is explicit and nbulk_bulk
    # is available; bulk and impurity modes have different semantics.
    if calculation_type in {"S", "L"}:
        if nlay is None or nlay <= 0:
            diagnostics.append(
                _error(
                    ErrorCode.COUNT_MISMATCH,
                    f"{calculation_type} mode requires a positive lattice nlay",
                    path=str(path),
                    field="nlay",
                    mode=calculation_type,
                    value=nlay,
                )
            )
        elif recursion_count is not None and recursion_count != nlay:
            diagnostics.append(
                _error(
                    ErrorCode.COUNT_MISMATCH,
                    f"{calculation_type} mode requires nrec=nlay={nlay}",
                    path=str(path),
                    field="nrec",
                    expected=nlay,
                    actual=recursion_count,
                    mode=calculation_type,
                )
            )
        if counts.get("nbulk_bulk") is not None and nlay is not None and input_ntype is not None:
            multiplier = 2 if calculation_type == "L" else 1
            expected_ntype = multiplier * counts["nbulk_bulk"] + nlay
            if input_ntype != expected_ntype:
                diagnostics.append(
                    _error(
                        ErrorCode.COUNT_MISMATCH,
                        f"{calculation_type} mode requires ntype={expected_ntype} from lattice counts",
                        path=str(path),
                        field="ntype",
                        expected=expected_ntype,
                        actual=input_ntype,
                        mode=calculation_type,
                    )
                )
    return normalized


def _validate_input_groups(
    lattice_group: Mapping[str, Any],
    atoms_group: Mapping[str, Any],
    lattice_path: Path,
    input_path: Path,
    lattice_data: Mapping[str, Any],
    diagnostics: list[Diagnostic],
) -> dict[str, Any]:
    normalized: dict[str, Any] = {}
    alat = _value(lattice_group, "alat")
    if alat is None:
        diagnostics.append(_error(ErrorCode.NAMELIST_FIELD_MISSING, "input &lattice field 'alat' is missing", path=str(input_path)))
    elif not _finite_number(alat) or float(alat) <= 0:
        diagnostics.append(_error(ErrorCode.NONPOSITIVE_SCALE, "input &lattice alat must be finite and positive", path=str(input_path), value=alat))
    else:
        normalized["alat"] = float(alat)

    labels = _labels(_value(atoms_group, "label"))
    if not labels:
        diagnostics.append(_error(ErrorCode.NAMELIST_FIELD_MISSING, "input &atoms must contain at least one label", path=str(input_path)))
    normalized["labels"] = labels

    explicit_ntype = _value(lattice_group, "ntype")
    ntype = _as_int(explicit_ntype) if explicit_ntype is not None else len(labels)
    if explicit_ntype is None:
        diagnostics.append(
            _info(
                ErrorCode.COUNT_MISMATCH,
                "input &lattice ntype is absent; derived from the label count for legacy compatibility",
                path=str(input_path),
                derived_ntype=ntype,
            )
        )
    elif ntype is None or ntype <= 0:
        diagnostics.append(_error(ErrorCode.COUNT_MISMATCH, "input &lattice ntype must be a positive integer", path=str(input_path)))
    if ntype is not None and ntype != len(labels):
        diagnostics.append(
            _error(
                ErrorCode.LABEL_COUNT_MISMATCH,
                f"input ntype={ntype} does not match {len(labels)} atom labels",
                path=str(input_path),
                expected=ntype,
                actual=len(labels),
            )
        )
    for index, label in enumerate(labels, start=1):
        symbol, atomic_number = _label_element(label)
        if symbol is None:
            diagnostics.append(
                _warning(
                    ErrorCode.LABEL_ELEMENT_UNKNOWN,
                    f"could not determine an element symbol from label '{label}'",
                    path=f"{input_path}:&atoms.label({index})",
                )
            )
        else:
            diagnostics.append(
                _info(
                    ErrorCode.LABEL_ELEMENT_DETERMINED,
                    f"label '{label}' resolved to element {symbol}",
                    path=f"{input_path}:&atoms.label({index})",
                    element=symbol,
                    atomic_number=atomic_number,
                )
            )
    normalized["ntype"] = ntype
    return normalized


def validate_setup_namelists(
    lattice_path: str | Path,
    input_path: str | Path,
    *,
    spec: SetupSpec | Mapping[str, Any] | None = None,
    generated: bool = False,
) -> ValidationResult:
    """Validate setup namelists and return coded semantic diagnostics."""

    lattice_file = Path(lattice_path)
    input_file = Path(input_path)
    diagnostics: list[Diagnostic] = []
    lattice_document = _read_namelist(lattice_file, diagnostics)
    input_document = _read_namelist(input_file, diagnostics)
    if lattice_document is None or input_document is None:
        return ValidationResult(valid=not _has_errors(diagnostics), diagnostics=tuple(diagnostics))

    lattice_group = _group(lattice_document, "lattice", diagnostics, lattice_file)
    input_lattice_group = _group(input_document, "lattice", diagnostics, input_file)
    atoms_group = _group(input_document, "atoms", diagnostics, input_file)
    if lattice_group is None or input_lattice_group is None or atoms_group is None:
        return ValidationResult(valid=not _has_errors(diagnostics), diagnostics=tuple(diagnostics))

    control_group = next(
        (
            value
            for key, value in input_document.items()
            if str(key).lower() == "control" and isinstance(value, Mapping)
        ),
        None,
    )
    calculation_type = None
    if control_group is not None:
        raw_calculation_type = _value(control_group, "calctype")
        if raw_calculation_type is not None:
            calculation_type = str(raw_calculation_type).strip().upper()

    raw_nlay = _value(input_lattice_group, "nlay")
    if raw_nlay is None:
        raw_nlay = _value(lattice_group, "nlay")
    nlay = _as_int(raw_nlay) if raw_nlay is not None else None
    raw_ntype = _value(input_lattice_group, "ntype")
    input_ntype = _as_int(raw_ntype) if raw_ntype is not None else len(_labels(_value(atoms_group, "label")))

    normalized_lattice = _validate_lattice_group(
        lattice_group,
        lattice_file,
        diagnostics,
        calculation_type=calculation_type,
        nlay=nlay,
        input_ntype=input_ntype,
    )
    normalized_input = _validate_input_groups(
        input_lattice_group,
        atoms_group,
        lattice_file,
        input_file,
        normalized_lattice,
        diagnostics,
    )

    setup_spec: SetupSpec | None
    if spec is None:
        setup_spec = None
    elif isinstance(spec, SetupSpec):
        setup_spec = spec
    else:
        try:
            setup_spec = SetupSpec.from_dict(spec)
        except Exception as exc:
            setup_spec = None
            diagnostics.append(_error(ErrorCode.SETUP_SPEC_INVALID, f"invalid setup specification: {exc}"))

    if setup_spec is not None:
        stochastic = _stochastic_provenance(setup_spec)
        if not stochastic["reproducible"]:
            diagnostics.append(
                _warning(
                    ErrorCode.UNSEEDED_SETUP,
                    f"setup contains {stochastic['reason']}; reproducibility is not guaranteed",
                )
            )
        if setup_spec.supercell != (1, 1, 1):
            expected_factor = math.prod(setup_spec.supercell)
            if expected_factor < 1:
                diagnostics.append(_error(ErrorCode.SUPERCELL_INVALID, "supercell dimensions must be positive"))
            normalized_input["supercell"] = list(setup_spec.supercell)

    if generated:
        control_group = _group(input_document, "control", diagnostics, input_file)
        if control_group is None:
            diagnostics.append(_error(ErrorCode.RECURSION_INVALID, "generated input must contain &control", path=str(input_file)))
        else:
            method = _value(control_group, "recur")
            steps = _as_int(_value(control_group, "lld"))
            if str(method).strip().lower() not in {"block", "chebyshev", "lanczos"}:
                diagnostics.append(_error(ErrorCode.RECURSION_INVALID, "control recur must be block, chebyshev, or lanczos", path=str(input_file)))
            if steps is None or steps <= 0:
                diagnostics.append(_error(ErrorCode.RECURSION_INVALID, "control lld must be a positive integer", path=str(input_file)))
            if setup_spec is not None:
                if str(method).strip().lower() != setup_spec.rec_method:
                    diagnostics.append(_warning(ErrorCode.RECURSION_INVALID, "generated recursion method differs from setup specification; preserving legacy output", path=str(input_file), expected=setup_spec.rec_method, actual=method))
                if steps != setup_spec.rec_steps:
                    diagnostics.append(_warning(ErrorCode.RECURSION_INVALID, "generated recursion steps differ from setup specification; preserving legacy output", path=str(input_file), expected=setup_spec.rec_steps, actual=steps))

        if setup_spec is not None:
            if setup_spec.preset == "KS":
                self_group = _group(input_document, "self", diagnostics, input_file)
                if self_group is not None and _value(self_group, "use_kspace") is not True:
                    diagnostics.append(
                        _error(
                            ErrorCode.SETUP_SPEC_INVALID,
                            "KS preset output must set &self use_kspace = .true.",
                            path=str(input_file),
                        )
                    )
            if setup_spec.mode_2d:
                b3 = _value(input_lattice_group, "b3")
                n3 = _as_int(_value(input_lattice_group, "n3"))
                if b3 is not False or n3 != 1:
                    diagnostics.append(_error(ErrorCode.TWO_D_INVALID, "2D setup output must disable b3 and set n3=1", path=str(input_file), b3=b3, n3=n3))
                vectors = normalized_lattice.get("a", [])
                if len(vectors) == 3 and _determinant(vectors) is not None and abs(float(_determinant(vectors))) <= 1.0e-10:
                    diagnostics.append(_error(ErrorCode.TWO_D_INVALID, "2D setup output has a singular cell", path=str(lattice_file)))

    return ValidationResult(valid=not _has_errors(diagnostics), diagnostics=tuple(diagnostics))


def validate_setup_outputs(
    output_dir: str | Path,
    *,
    lattice_path: str | Path | None = None,
    input_path: str | Path | None = None,
    spec: SetupSpec | Mapping[str, Any] | None = None,
    check_supercell_expansion: bool = False,
) -> ValidationResult:
    """Validate generated setup files, including sphere/radius evidence."""

    directory = Path(output_dir)
    diagnostics: list[Diagnostic] = []
    generated_lattice = directory / "lattice_new.nml"
    generated_input = directory / "input_new.nml"
    radii_path = directory / "radii.out"
    for path in (generated_lattice, generated_input, radii_path):
        if not path.is_file():
            diagnostics.append(_error(ErrorCode.OUTPUT_MISSING, f"expected generated output is missing: {path.name}", path=str(path)))
        elif not path.stat().st_size:
            diagnostics.append(_error(ErrorCode.OUTPUT_UNREADABLE, f"generated output is empty: {path.name}", path=str(path)))

    if generated_lattice.is_file() and generated_input.is_file():
        pair_validation = validate_setup_namelists(generated_lattice, generated_input, spec=spec, generated=True)
        diagnostics.extend(pair_validation.diagnostics)

    if radii_path.is_file() and radii_path.stat().st_size:
        _validate_radii(radii_path, generated_lattice, diagnostics)

    if lattice_path is not None and input_path is not None:
        source_validation = validate_setup_namelists(lattice_path, input_path, spec=spec)
        diagnostics.extend(source_validation.diagnostics)
        if check_supercell_expansion and isinstance(spec, SetupSpec) and spec.supercell != (1, 1, 1):
            _validate_supercell_expansion(
                Path(lattice_path),
                Path(input_path),
                generated_lattice,
                generated_input,
                spec,
                diagnostics,
            )
    return ValidationResult(valid=not _has_errors(diagnostics), diagnostics=tuple(diagnostics))


def _validate_supercell_expansion(
    source_lattice: Path,
    source_input: Path,
    generated_lattice: Path,
    generated_input: Path,
    spec: SetupSpec,
    diagnostics: list[Diagnostic],
) -> None:
    """Check the count expansion without imposing species-order assumptions."""

    try:
        source_lattice_group = f90nml.read(source_lattice)["lattice"]
        generated_lattice_group = f90nml.read(generated_lattice)["lattice"]
        source_input_document = f90nml.read(source_input)
        generated_input_document = f90nml.read(generated_input)
        source_counts = {
            name: _as_int(_value(source_lattice_group, name))
            for name in ("nbulk_bulk", "ntot", "nbas", "nrec")
        }
        generated_counts = {
            name: _as_int(_value(generated_lattice_group, name))
            for name in ("nbulk_bulk", "ntot", "nbas", "nrec")
        }
        source_labels = _labels(_value(source_input_document["atoms"], "label"))
        generated_labels = _labels(_value(generated_input_document["atoms"], "label"))
    except Exception as exc:
        diagnostics.append(_error(ErrorCode.SUPERCELL_INVALID, f"could not validate supercell expansion: {exc}"))
        return
    source_ntot = source_counts.get("ntot")
    factor = math.prod(spec.supercell)
    for name, source_value in source_counts.items():
        generated_value = generated_counts.get(name)
        if source_value is None or generated_value is None:
            continue
        expected = source_value * factor if name == "ntot" or source_value == source_ntot else source_value
        if generated_value != expected:
            diagnostics.append(
                _error(
                    ErrorCode.SUPERCELL_INVALID,
                    f"generated {name} does not preserve its source-defined scaling",
                    field=name,
                    expected=expected,
                    actual=generated_value,
                )
            )
    if source_ntot is not None and len(generated_labels) != len(source_labels) * factor:
        diagnostics.append(_error(ErrorCode.SUPERCELL_INVALID, "generated label count does not match supercell factor", expected=len(source_labels) * math.prod(spec.supercell), actual=len(generated_labels)))


def _validate_radii(path: Path, lattice_path: Path, diagnostics: list[Diagnostic]) -> None:
    try:
        lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    except Exception as exc:
        diagnostics.append(_error(ErrorCode.OUTPUT_UNREADABLE, f"could not read radii output: {exc}", path=str(path)))
        return
    if not lines:
        diagnostics.append(_error(ErrorCode.SPHERE_OUTPUT_INVALID, "radii.out contains no records", path=str(path)))
        return
    count = _as_int(lines[0].split()[0]) if lines[0].split() else None
    records = lines[1:]
    if count is None or count <= 0:
        diagnostics.append(_error(ErrorCode.SPHERE_OUTPUT_INVALID, "radii.out must start with a positive atom count", path=str(path)))
        return
    if lattice_path.is_file():
        try:
            generated_group = f90nml.read(lattice_path)["lattice"]
            generated_count = _as_int(_value(generated_group, "nrec"))
            if generated_count is not None and count != generated_count:
                diagnostics.append(
                    _error(
                        ErrorCode.SPHERE_OUTPUT_INVALID,
                        "radii.out atom count does not match generated lattice nrec",
                        path=str(path),
                        expected=generated_count,
                        actual=count,
                    )
                )
        except Exception:
            pass
    if len(records) != count:
        diagnostics.append(_error(ErrorCode.SPHERE_OUTPUT_INVALID, f"radii.out declares {count} records but contains {len(records)}", path=str(path), expected=count, actual=len(records)))
    for index, line in enumerate(records, start=1):
        fields = line.replace("D", "E").replace("d", "e").split()
        if len(fields) < 5:
            diagnostics.append(_error(ErrorCode.SPHERE_OUTPUT_INVALID, f"radii record {index} has too few fields", path=str(path)))
            continue
        try:
            atom_index = int(fields[0])
            coordinates = [float(value) for value in fields[1:4]]
            radius = float(fields[4])
        except ValueError:
            diagnostics.append(_error(ErrorCode.SPHERE_OUTPUT_INVALID, f"radii record {index} is not numeric", path=str(path)))
            continue
        if atom_index <= 0 or not all(math.isfinite(value) for value in coordinates) or not math.isfinite(radius) or radius <= 0:
            diagnostics.append(_error(ErrorCode.SPHERE_OUTPUT_INVALID, f"radii record {index} contains invalid index, coordinate, or radius", path=str(path)))


def _now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _git_sha() -> str | None:
    try:
        root = Path(__file__).resolve().parent.parent
        result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, check=False)
        value = result.stdout.strip()
        return value or None
    except OSError:
        return None


def _fingerprint(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.is_file():
        return None
    try:
        return {"path": path.name, "sha256": hash_file(path), "size": path.stat().st_size}
    except OSError:
        return {"path": path.name}


def _software_metadata() -> dict[str, Any]:
    metadata: dict[str, Any] = {
        "rsutil_version": None,
        "git_sha": _git_sha(),
        "python_executable": _fingerprint(Path(sys.executable)),
    }
    try:
        import rsutil

        metadata["rsutil_version"] = rsutil.__version__
        library = getattr(rsutil, "lib_path", None)
        metadata["shared_library"] = _fingerprint(Path(library)) if library else None
    except Exception:
        metadata["shared_library"] = None
    return metadata


def _environment_metadata() -> dict[str, Any]:
    return {
        "platform": platform.platform(),
        "system": platform.system(),
        "machine": platform.machine(),
        "python": sys.version,
        "python_implementation": platform.python_implementation(),
        "python_compiler": platform.python_compiler(),
        "compiler": sysconfig.get_config_var("CC"),
    }


def _structure_identity(lattice_path: Path, input_path: Path) -> dict[str, Any] | None:
    try:
        lattice = f90nml.read(lattice_path)["lattice"]
        document = f90nml.read(input_path)
        input_lattice = document["lattice"]
        atoms = document["atoms"]
        return {
            "lattice": {
                "a": _value(lattice, "a"),
                "crd": _value(lattice, "crd"),
                "nrec": _value(lattice, "nrec"),
            },
            "input_lattice": {
                "alat": _value(input_lattice, "alat"),
                "pbc": _value(input_lattice, "pbc"),
                "b1": _value(input_lattice, "b1"),
                "b2": _value(input_lattice, "b2"),
                "b3": _value(input_lattice, "b3"),
            },
            "labels": _labels(_value(atoms, "label")),
        }
    except Exception:
        return None


def build_run_manifest(
    lattice_path: str | Path,
    input_path: str | Path,
    spec: SetupSpec | Mapping[str, Any] | None,
    *,
    status: OperationStatus | str,
    validation: ValidationResult | None = None,
    output_paths: Mapping[str, str | Path] | None = None,
    started_at: str | None = None,
    finished_at: str | None = None,
    base_provenance: Mapping[str, Any] | None = None,
) -> RunManifest:
    """Build a deterministic-identity manifest for setup inputs and outputs."""

    lattice_file = Path(lattice_path)
    input_file = Path(input_path)
    setup_spec = spec if isinstance(spec, SetupSpec) else SetupSpec.from_dict(spec or {})
    structure = _structure_identity(lattice_file, input_file)
    structure_hash = scientific_hash(structure) if structure is not None else None
    normalized_spec = setup_spec.to_dict()
    setup_identity_hash = scientific_hash({"structure_hash": structure_hash, "setup_spec": normalized_spec})
    input_hashes: dict[str, Any] = {}
    for name, path in (("lattice.nml", lattice_file), ("input.nml", input_file)):
        if path.is_file():
            input_hashes[name] = {"sha256": hash_file(path), "size": path.stat().st_size}
    output_hashes: dict[str, Any] = {}
    for name, path_value in (output_paths or {}).items():
        path = Path(path_value)
        if path.is_file():
            output_hashes[name] = {"sha256": hash_file(path), "size": path.stat().st_size}
    diagnostics = () if validation is None else validation.diagnostics
    validation_data = {
        "valid": validation.valid if validation is not None else status == OperationStatus.SUCCESS,
        "error_count": sum(item.severity is DiagnosticSeverity.ERROR for item in diagnostics),
        "warning_count": sum(item.severity is DiagnosticSeverity.WARNING for item in diagnostics),
        "diagnostics": [item.to_dict() for item in diagnostics],
    }
    stochastic = _stochastic_provenance(setup_spec)
    timestamps = {"started_at": started_at or _now(), "finished_at": finished_at or _now()}
    provenance = dict(base_provenance or {})
    provenance.pop("timestamps", None)
    software = _software_metadata()
    if provenance:
        software["worker"] = provenance
    return RunManifest(
        status=status,
        input_hashes=input_hashes,
        output_hashes=output_hashes,
        normalized_setup_spec=normalized_spec,
        structure_hash=structure_hash,
        setup_identity_hash=setup_identity_hash,
        software=software,
        environment=_environment_metadata(),
        timestamps=timestamps,
        validation=validation_data,
        stochastic=stochastic,
    )


def setup_report(result: SetupResult, manifest: RunManifest | None = None) -> SetupReport:
    """Wrap a setup result in the versioned report envelope."""

    if manifest is None:
        raw = result.provenance.get("manifest")
        if isinstance(raw, RunManifest):
            manifest = raw
        elif isinstance(raw, Mapping):
            manifest = RunManifest.from_dict(raw)
        else:
            manifest = RunManifest(status=result.status)
    return SetupReport(result=result, manifest=manifest)


def write_setup_report(report: SetupReport | SetupResult, path: str | Path) -> Path:
    """Serialize a report deterministically enough for machine consumption."""

    envelope = setup_report(report) if isinstance(report, SetupResult) else report
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(envelope.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return target


def read_setup_report(path: str | Path) -> SetupReport:
    """Read a versioned report into the typed report model."""

    with Path(path).open("r", encoding="utf-8") as stream:
        return SetupReport.from_dict(json.load(stream))


validate_setup = validate_setup_namelists


__all__ = [
    "build_run_manifest",
    "read_setup_report",
    "setup_report",
    "validate_setup",
    "validate_setup_namelists",
    "validate_setup_outputs",
    "write_setup_report",
]
