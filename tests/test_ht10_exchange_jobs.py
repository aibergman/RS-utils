"""HT-10 transactional exchange writer and manifest coverage."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from rsutil.exchange.jobs import (
    EXCHANGE_MANIFEST_FILENAME,
    EXCHANGE_MANIFEST_SCHEMA,
    ExchangeCollisionError,
    ExchangeManifestError,
    read_exchange_manifest,
)
from rsutil.exchange.service import ExchangeService
from rsutil.core.hashing import hash_file

from .conftest import copy_example_files


def _copy_exchange_fixture(destination: Path) -> None:
    copy_example_files(
        destination,
        "FePt",
        [
            "input.nml",
            "lattice.nml",
            "clust",
            "report.out",
            "Fe.nml",
            "Pt.nml",
            "Fe_out.nml",
            "Pt_out.nml",
        ],
    )


def test_manifest_is_canonical_and_identical_reuse_is_non_mutating(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)
    service = ExchangeService(tmp_path)

    first = service.generate(cutoff=1.0)
    manifest_path = tmp_path / "Atom_1" / EXCHANGE_MANIFEST_FILENAME
    before = manifest_path.read_bytes()
    manifest = read_exchange_manifest(manifest_path)

    assert manifest["schema_version"] == EXCHANGE_MANIFEST_SCHEMA
    assert manifest["canonical_pair_list"] == sorted(
        manifest["canonical_pair_list"],
        key=lambda item: (item["i"], item["j"], *item["R"]),
    )
    assert manifest["spec_hash"]
    assert manifest["source_run"]["hash"] == first.jobs[0].source_run_hash

    second = service.generate(cutoff=1.0, collision_policy="reuse_if_identical")
    assert second.write_result is not None
    assert tmp_path / "Atom_1" in second.write_result.reused
    assert manifest_path.read_bytes() == before


def test_runtime_outputs_are_ignored_during_identical_reuse(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)
    service = ExchangeService(tmp_path)
    service.generate(cutoff=1.0)
    runtime_output = tmp_path / "Atom_1" / "runtime.log"
    runtime_output.write_text("solver log\n", encoding="utf-8")

    result = service.generate(cutoff=1.0, collision_policy="reuse_if_identical")

    assert result.write_result is not None
    assert tmp_path / "Atom_1" in result.write_result.reused
    assert runtime_output.read_text(encoding="utf-8") == "solver log\n"


def test_reuse_rejects_coordinated_input_and_historical_manifest_tampering(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)
    service = ExchangeService(tmp_path)
    service.generate(cutoff=1.0)
    atom = tmp_path / "Atom_1"
    manifest_path = atom / EXCHANGE_MANIFEST_FILENAME
    tampered = atom / "input.nml"
    tampered.write_text(tampered.read_text(encoding="utf-8") + "\n! tampered\n", encoding="utf-8")

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["files"]["input.nml"] = hash_file(tampered)
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ExchangeCollisionError, match="freshly generated inputs"):
        service.generate(cutoff=1.0, collision_policy="reuse_if_identical")


def test_changed_pair_set_rejects_stale_reuse_and_does_not_replace_job(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)
    service = ExchangeService(tmp_path)
    service.generate(cutoff=1.0)
    manifest_path = tmp_path / "Atom_1" / EXCHANGE_MANIFEST_FILENAME
    before = manifest_path.read_bytes()

    with pytest.raises(ExchangeCollisionError, match="stale reuse|not identical"):
        service.generate(cutoff=0.5, collision_policy="reuse_if_identical")

    assert manifest_path.read_bytes() == before


def test_unmarked_atom_directory_is_never_deleted_by_default(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)
    atom = tmp_path / "Atom_1"
    atom.mkdir()
    user_file = atom / "user_data.txt"
    user_file.write_text("keep me", encoding="utf-8")

    with pytest.raises(ExchangeCollisionError, match="unmarked"):
        ExchangeService(tmp_path).generate(cutoff=1.0, collision_policy="replace_generated")

    assert user_file.read_text(encoding="utf-8") == "keep me"
    assert not (atom / EXCHANGE_MANIFEST_FILENAME).exists()


def test_staging_failure_preserves_existing_generated_job(tmp_path: Path, monkeypatch) -> None:
    _copy_exchange_fixture(tmp_path)
    service = ExchangeService(tmp_path)
    service.generate(cutoff=1.0)
    atom = tmp_path / "Atom_1"
    before = {
        path.name: path.read_bytes()
        for path in atom.iterdir()
        if path.is_file()
    }

    import rsutil.exchange.jobs as jobs_module

    def fail_write(*args, **kwargs):
        raise OSError("simulated staged write failure")

    monkeypatch.setattr(jobs_module.f90nml, "write", fail_write)
    with pytest.raises(OSError, match="simulated"):
        service.generate(cutoff=1.0, collision_policy="replace_generated")

    assert {
        path.name: path.read_bytes()
        for path in atom.iterdir()
        if path.is_file()
    } == before


def test_manifest_ijpair_matches_generated_input(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)
    ExchangeService(tmp_path).generate(cutoff=1.0)

    for atom_index in (1, 2):
        atom = tmp_path / f"Atom_{atom_index}"
        manifest = json.loads(
            (atom / EXCHANGE_MANIFEST_FILENAME).read_text(encoding="utf-8")
        )
        import f90nml

        generated = f90nml.read(atom / "input.nml")
        assert manifest["legacy_ijpair"] == generated["lattice"]["ijpair"]
        assert generated["lattice"]["njij"] == len(manifest["canonical_pair_list"])


def test_tampered_compatibility_pair_view_is_rejected(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)
    ExchangeService(tmp_path).generate(cutoff=1.0)
    manifest_path = tmp_path / "Atom_1" / EXCHANGE_MANIFEST_FILENAME
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["canonical_pair_list"][0]["R"][0] += 1
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ExchangeManifestError, match="canonical_pair_list|pair-list"):
        read_exchange_manifest(manifest_path)


def test_tampered_authoritative_pair_list_without_hash_update_is_rejected(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)
    ExchangeService(tmp_path).generate(cutoff=1.0)
    manifest_path = tmp_path / "Atom_1" / EXCHANGE_MANIFEST_FILENAME
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["spec"]["pair_keys"][0]["R"][0] += 1
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ExchangeManifestError, match="specification hash"):
        read_exchange_manifest(manifest_path)
