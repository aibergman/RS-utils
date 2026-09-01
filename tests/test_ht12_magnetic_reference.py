"""HT-R12 source-derived magnetic reference and symmetry contracts."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from rsutil.exchange.service import ExchangeService
from rsutil.exchange.symmetry import (
    MAGNETIC_REFERENCE_COLLINEAR_FM,
    MAGNETIC_REFERENCE_UNKNOWN,
)

from .test_ht08_exchange_service import _copy_exchange_fixture


def _reference_with_nsp(destination: Path, nsp: int | None) -> None:
    _copy_exchange_fixture(destination)
    if nsp is not None:
        input_path = destination / "input.nml"
        contents = input_path.read_text(encoding="utf-8")
        contents = re.sub(r"(?im)^(\s*nsp\s*=\s*)\d+", rf"\g<1>{nsp}", contents, count=1)
        input_path.write_text(contents, encoding="utf-8")


@pytest.mark.parametrize(
    ("nsp", "is_collinear", "is_noncollinear", "has_soc"),
    [
        (1, True, False, False),
        (2, True, False, True),
        (3, False, True, False),
        (4, False, True, True),
    ],
)
def test_reference_derives_nsp_spin_and_soc_contract(
    tmp_path: Path,
    nsp: int,
    is_collinear: bool,
    is_noncollinear: bool,
    has_soc: bool,
) -> None:
    _reference_with_nsp(tmp_path, nsp)

    reference = ExchangeService(tmp_path).plan().reference

    assert reference.nsp == nsp
    assert reference.is_collinear is is_collinear
    assert reference.is_noncollinear is is_noncollinear
    assert reference.has_soc is has_soc
    assert reference.nsp_semantics.known


def test_missing_nsp_is_explicitly_unknown(tmp_path: Path) -> None:
    _reference_with_nsp(tmp_path, None)
    input_path = tmp_path / "input.nml"
    contents = input_path.read_text(encoding="utf-8")
    contents = re.sub(r"(?im)^\s*nsp\s*=.*$\n?", "", contents, count=1)
    input_path.write_text(contents, encoding="utf-8")

    reference = ExchangeService(tmp_path).plan().reference

    assert reference.nsp is None
    assert reference.is_collinear is None
    assert reference.is_noncollinear is None
    assert reference.has_soc is None


def test_nsp_does_not_infer_ferromagnetic_order_or_reduce(tmp_path: Path) -> None:
    _reference_with_nsp(tmp_path, 1)

    with pytest.warns(UserWarning, match="identity-only"):
        unknown = ExchangeService(tmp_path).plan(cutoff=1.0, unique=True)
    full = ExchangeService(tmp_path).plan(cutoff=1.0, full_pairs=True)

    assert len(unknown.pairs) == len(full.pairs)
    assert all(job.magnetic_reference == MAGNETIC_REFERENCE_UNKNOWN for job in unknown.jobs)
    assert all(job.soc is False for job in unknown.jobs)
    assert all(job.symmetry_reduction_applied is False for job in unknown.jobs)
    assert all(job.symmetry_reduction_reason for job in unknown.jobs)


def test_nsp_two_enables_soc_and_rejects_conflicting_explicit_default(
    tmp_path: Path,
) -> None:
    _reference_with_nsp(tmp_path, 2)
    service = ExchangeService(tmp_path)

    plan = service.plan()
    assert plan.reference.has_soc is True
    assert all(job.soc is True for job in plan.jobs)
    assert all(job.soc_source == "reference:control.nsp" for job in plan.jobs)

    with pytest.raises(ValueError, match="conflicts with reference nsp=2"):
        service.plan(soc=False)


@pytest.mark.parametrize("nsp", [3, 4])
def test_noncollinear_nsp_is_identity_only_and_never_fm(
    tmp_path: Path, nsp: int
) -> None:
    _reference_with_nsp(tmp_path, nsp)

    with pytest.warns(UserWarning, match="identity-only"):
        plan = ExchangeService(tmp_path).plan(cutoff=1.0, unique=True)

    assert all(job.nsp == nsp for job in plan.jobs)
    assert all(job.is_noncollinear is True for job in plan.jobs)
    assert all(job.magnetic_reference == MAGNETIC_REFERENCE_UNKNOWN for job in plan.jobs)
    assert all(job.soc is (nsp == 4) for job in plan.jobs)
    assert all(job.symmetry_reduction_applied is False for job in plan.jobs)


def test_explicit_fm_no_soc_enables_only_proven_reduction(tmp_path: Path) -> None:
    _reference_with_nsp(tmp_path, 1)
    service = ExchangeService(tmp_path)

    reduced = service.plan(
        cutoff=1.0,
        unique=True,
        magnetic_reference=MAGNETIC_REFERENCE_COLLINEAR_FM,
    )
    full = service.plan(cutoff=1.0, full_pairs=True)

    assert len(reduced.pairs) < len(full.pairs)
    assert all(job.magnetic_reference == MAGNETIC_REFERENCE_COLLINEAR_FM for job in reduced.jobs)
    assert all(job.symmetry_reduction_applied is True for job in reduced.jobs)


def test_explicit_fm_soc_requires_axis_and_valid_axis_is_restricted(
    tmp_path: Path,
) -> None:
    _reference_with_nsp(tmp_path, 2)
    service = ExchangeService(tmp_path)
    full = service.plan(cutoff=1.0, full_pairs=True)

    with pytest.warns(UserWarning, match="explicit global magnetization axis"):
        no_axis = service.plan(
            cutoff=1.0,
            unique=True,
            magnetic_reference=MAGNETIC_REFERENCE_COLLINEAR_FM,
        )
    assert len(no_axis.pairs) == len(full.pairs)
    assert all(job.symmetry_reduction_applied is False for job in no_axis.jobs)

    with_axis = service.plan(
        cutoff=1.0,
        unique=True,
        magnetic_reference=MAGNETIC_REFERENCE_COLLINEAR_FM,
        magnetic_axis=(0.0, 0.0, 1.0),
    )
    assert len(with_axis.pairs) < len(full.pairs)
    assert all(job.symmetry_reduction_applied is True for job in with_axis.jobs)


def test_nsp_four_cannot_be_overridden_to_collinear_fm(tmp_path: Path) -> None:
    _reference_with_nsp(tmp_path, 4)

    with pytest.raises(ValueError, match="conflicts with noncollinear nsp=4"):
        ExchangeService(tmp_path).plan(
            unique=True,
            magnetic_reference=MAGNETIC_REFERENCE_COLLINEAR_FM,
            magnetic_axis=(0.0, 0.0, 1.0),
        )


def test_exchange_manifest_records_source_magnetic_context(tmp_path: Path) -> None:
    _reference_with_nsp(tmp_path, 4)

    result = ExchangeService(tmp_path).generate(cutoff=1.0, unique=True)
    manifest_path = result.jobs[0].output_dir / "rsutil_exchange_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    spec = manifest["spec"]

    assert spec["nsp"] == 4
    assert spec["nsp_source"] == "input.nml"
    assert spec["is_collinear"] is False
    assert spec["is_noncollinear"] is True
    assert spec["soc_source"] == "reference:control.nsp"
    assert spec["magnetic_reference"] == MAGNETIC_REFERENCE_UNKNOWN
    assert spec["symmetry_reduction_applied"] is False
    assert "noncollinear" in spec["symmetry_reduction_reason"]
