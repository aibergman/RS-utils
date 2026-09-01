"""Focused HT-02 contract and safe-path tests."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from rsutil.core import (
    CoordinateConvention,
    RSStructure,
    SetupSpec,
    StructureSite,
    ValidationError,
    Workspace,
    atomic_replace,
    canonical_json,
    hash_file,
    isolated_workspace,
    scientific_hash,
    staging_workspace,
)


def _structure() -> RSStructure:
    return RSStructure(
        lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        length_unit="angstrom",
        coordinate_convention=CoordinateConvention.FRACTIONAL,
        sites=(
            StructureSite(0, (0, 0, 0), "Fe", 26, rs_label="Fe1", rs_type="host"),
            StructureSite(1, (0.5, 0.5, 0.5), "N", 7, rs_label="N1", rs_type="interstitial"),
        ),
    )


def test_structure_validation_and_round_trip() -> None:
    structure = _structure()
    restored = RSStructure.from_dict(structure.to_dict())
    assert restored == structure
    assert restored.positions == ((0.0, 0.0, 0.0), (0.5, 0.5, 0.5))
    with pytest.raises(ValidationError):
        RSStructure(
            lattice_vectors=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            length_unit="angstrom",
            sites=(StructureSite(1, (0, 0, 0), "Fe", 26),),
        )


def test_setup_spec_is_scientific_only_and_validates() -> None:
    spec = SetupSpec(supercell=(2, 1, 1), seed=42, sqs="Fe=0.5,Co=0.5")
    assert spec.to_dict()["supercell"] == [2, 1, 1]
    assert not hasattr(spec, "in_place")
    with pytest.raises(ValidationError):
        SetupSpec(z_iso_factor=2.0)


def test_hashing_is_order_independent_and_excludes_timestamps() -> None:
    left = {"b": [2, 3], "a": {"value": 1, "timestamp": "old"}}
    right = {"a": {"timestamp": "new", "value": 1}, "b": [2, 3]}
    assert canonical_json(left) != canonical_json(right)
    assert scientific_hash(left) == scientific_hash(right)


def test_file_hash_is_content_only(tmp_path: Path) -> None:
    source = tmp_path / "source.dat"
    source.write_bytes(b"abc")
    assert hash_file(source) == hashlib.sha256(b"abc").hexdigest()


def test_workspace_cleans_up_after_exception(tmp_path: Path) -> None:
    holder: dict[str, Path] = {}
    with pytest.raises(RuntimeError):
        with isolated_workspace(tmp_path) as workspace:
            holder["path"] = workspace
            (workspace / "partial").write_text("candidate", encoding="utf-8")
            raise RuntimeError("failure")
    assert not holder["path"].exists()


def test_workspace_child_rejects_absolute_dotdot_and_symlink_escape(tmp_path: Path) -> None:
    workspace = Workspace.create(tmp_path)
    outside = tmp_path / "outside"
    outside.mkdir()
    link = workspace.path / "link"
    link.symlink_to(outside, target_is_directory=True)
    try:
        for name in ("../escape", str(tmp_path / "absolute"), "link/escape"):
            with pytest.raises(ValueError, match="escapes workspace"):
                workspace.child(name)
        assert workspace.child("nested/file").is_relative_to(workspace.path.resolve())
    finally:
        workspace.cleanup()


def test_workspaces_are_independent_and_same_filesystem_staging_is_atomic(tmp_path: Path) -> None:
    first = Workspace.create(tmp_path)
    second = Workspace.create(tmp_path)
    try:
        assert first.path != second.path
        (first.path / "marker").write_text("first", encoding="utf-8")
        assert not (second.path / "marker").exists()
    finally:
        first.cleanup()
        second.cleanup()

    destination = tmp_path / "result.txt"
    destination.write_text("old", encoding="utf-8")
    with staging_workspace(destination) as staging:
        candidate = staging / "candidate.txt"
        candidate.write_text("new", encoding="utf-8")
        atomic_replace(candidate, destination)
        assert destination.read_text(encoding="utf-8") == "new"
    assert not staging.exists()
