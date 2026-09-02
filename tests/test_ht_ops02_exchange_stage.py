"""HT-OPS-02 campaign exchange-stage coverage."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

from rsutil.campaign import (
    EXECUTION_EXITED_ZERO,
    EXCHANGE_STAGE_MANIFEST_FILENAME,
    ExchangeStagePreparation,
    CampaignIndex,
    CampaignIndexEntry,
    COMPLETED_INVALID,
    OUTPUT_NO_OUTPUT,
    collect_exchange_stage,
    prepare_exchange_stage,
    render_exchange_slurm,
    status_exchange_stage,
)

from .conftest import copy_example_files


def _write_parent(
    root: Path,
    *,
    entry_id: str,
    name: str,
    report: str = "Fermi energy: -0.125\n",
    execution: bool = True,
) -> CampaignIndexEntry:
    entry = CampaignIndexEntry(
        entry_id,
        name,
        f"runs/{name}",
        f"runs/{name}/run_manifest.json",
        f"manifest-{name}",
    )
    run = root / entry.run_path
    run.mkdir(parents=True)
    copy_example_files(
        run,
        "FePt",
        [
            "lattice.nml",
            "input.nml",
            "clust",
            "Fe_out.nml",
            "Pt_out.nml",
        ],
    )
    (run / "report.out").write_text(report, encoding="utf-8")
    (run / "out.scf").write_text("Converged!\n" if "Not converged" not in report else "Not converged!\n", encoding="utf-8")
    (run / "run_manifest.json").write_text(
        json.dumps({
            "schema_version": "rsutil-run-v2",
            "campaign": {
                "entry_id": entry.entry_id,
                "manifest_hash": entry.manifest_hash,
            },
        }),
        encoding="utf-8",
    )
    if execution:
        (run / "rsutil_execution.json").write_text(
            json.dumps({
                "schema_version": "rsutil-execution-v1",
                "identity": {
                    "entry_id": entry.entry_id,
                    "manifest_hash": entry.manifest_hash,
                },
                "status": "success",
                "return_code": 0,
            }),
            encoding="utf-8",
        )
    return entry


def _write_campaign(root: Path, entries: tuple[CampaignIndexEntry, ...]) -> None:
    index = CampaignIndex("campaign-id", "exchange-parent", "spec-hash", entries)
    (root / "campaign_index.json").write_text(
        json.dumps(index.to_dict()), encoding="utf-8"
    )


def test_prepare_selects_only_converged_successful_parents_and_reuses_jobs(tmp_path: Path) -> None:
    valid = _write_parent(tmp_path, entry_id="a" * 64, name="valid")
    nonconverged = _write_parent(
        tmp_path,
        entry_id="b" * 64,
        name="nonconverged",
        report="Fermi energy: -0.125\nNot converged!\n",
    )
    failed = _write_parent(
        tmp_path,
        entry_id="c" * 64,
        name="failed",
        execution=False,
    )
    _write_campaign(tmp_path, (valid, nonconverged, failed))

    first = prepare_exchange_stage(tmp_path, tmp_path / "exchange-stage", cutoff=1.0)
    assert isinstance(first, ExchangeStagePreparation)
    payload = json.loads((tmp_path / "exchange-stage" / EXCHANGE_STAGE_MANIFEST_FILENAME).read_text())
    assert payload["parent_entry_ids"] == [valid.entry_id]
    assert payload["parent_run_paths"] == [valid.run_path]
    assert payload["exchange_options"]["pair_policy"] == "full"
    assert len(first.jobs) == 2
    assert {item["name"] for item in first.skipped_parents} == {"nonconverged", "failed"}

    before = {
        item["exchange_manifest_path"]: (tmp_path / "exchange-stage" / item["exchange_manifest_path"]).read_bytes()
        for item in first.jobs
    }
    second = prepare_exchange_stage(tmp_path, tmp_path / "exchange-stage", cutoff=1.0)
    assert second.stage_id == first.stage_id
    assert {
        item["exchange_manifest_path"]: (tmp_path / "exchange-stage" / item["exchange_manifest_path"]).read_bytes()
        for item in second.jobs
    } == before

    changed = prepare_exchange_stage(tmp_path, tmp_path / "exchange-stage", cutoff=0.5)
    assert changed.stage_id != first.stage_id
    assert changed.jobs[0]["canonical_exchange_job_identity"] != first.jobs[0]["canonical_exchange_job_identity"]


def test_exchange_stage_identity_tracks_parent_output_and_persists_each_stage(tmp_path: Path) -> None:
    entry = _write_parent(tmp_path, entry_id="a" * 64, name="valid")
    _write_campaign(tmp_path, (entry,))

    first = prepare_exchange_stage(tmp_path, tmp_path / "exchange-stage", cutoff=1.0)
    first_manifest = first.manifest_path.read_bytes()

    # The parent entry and setup manifest are unchanged, but the converged
    # exchange source evidence is different.
    (tmp_path / entry.run_path / "report.out").write_text(
        "Fermi energy: -0.250\n", encoding="utf-8"
    )
    changed_parent = prepare_exchange_stage(
        tmp_path, tmp_path / "exchange-stage", cutoff=1.0
    )
    assert changed_parent.stage_id != first.stage_id
    assert first.manifest_path.is_file()
    assert first.manifest_path.read_bytes() == first_manifest
    assert changed_parent.manifest_path.is_file()

    changed_policy = prepare_exchange_stage(
        tmp_path, tmp_path / "exchange-stage", cutoff=0.5
    )
    assert changed_policy.stage_id != changed_parent.stage_id
    assert changed_policy.manifest_path.is_file()
    assert len(list((tmp_path / "exchange-stage").glob("*/exchange_stage_manifest.json"))) == 3


def test_unique_reduction_requires_explicit_supported_magnetic_context(tmp_path: Path) -> None:
    entry = _write_parent(tmp_path, entry_id="a" * 64, name="valid")
    _write_campaign(tmp_path, (entry,))

    try:
        prepare_exchange_stage(
            tmp_path,
            tmp_path / "exchange-stage",
            unique=True,
            full_pairs=False,
        )
    except ValueError as exc:
        assert "magnetic_reference" in str(exc)
    else:  # pragma: no cover - assertion is the meaningful branch
        raise AssertionError("unsupported magnetic context was accepted")


def test_identical_parent_runs_keep_distinct_stable_array_identities(tmp_path: Path) -> None:
    first = _write_parent(tmp_path, entry_id="a" * 64, name="first")
    second = _write_parent(tmp_path, entry_id="b" * 64, name="second")
    _write_campaign(tmp_path, (first, second))

    stage = prepare_exchange_stage(tmp_path, tmp_path / "exchange-stage")
    assert len(stage.jobs) == 4
    assert len({item["job_identity"] for item in stage.jobs}) == 4
    assert {item["parent_entry_id"] for item in stage.jobs} == {
        first.entry_id,
        second.entry_id,
    }


def test_explicit_collinear_fm_soc_policy_enables_unique_reduction(tmp_path: Path) -> None:
    entry = _write_parent(tmp_path, entry_id="a" * 64, name="valid")
    _write_campaign(tmp_path, (entry,))

    stage = prepare_exchange_stage(
        tmp_path,
        tmp_path / "exchange-stage",
        unique=True,
        full_pairs=False,
        soc=True,
        magnetic_reference="collinear_fm",
        magnetic_axis=(0.0, 0.0, 1.0),
    )
    assert stage.jobs
    payload = json.loads((stage.manifest_path).read_text())
    assert payload["exchange_options"]["pair_policy"] == "unique"
    assert all(
        json.loads((stage.stage_root / item["exchange_manifest_path"]).read_text())["spec"]["symmetry_reduction_applied"]
        for item in stage.jobs
    )


def test_exchange_array_records_execution_evidence_and_keeps_missing_output_invalid(tmp_path: Path) -> None:
    entry = _write_parent(tmp_path, entry_id="a" * 64, name="valid")
    _write_campaign(tmp_path, (entry,))
    stage = prepare_exchange_stage(tmp_path, tmp_path / "exchange-stage")
    fake = tmp_path / "fake exchange.sh"
    fake.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    fake.chmod(0o755)
    rendered = render_exchange_slurm(stage.stage_root, executable=fake)
    environment = os.environ.copy()
    environment["SLURM_ARRAY_TASK_ID"] = "0"
    completed = subprocess.run(
        ["bash", str(rendered.script_path)],
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    report = status_exchange_stage(stage.stage_root)
    assert report.entries[0].execution_status == EXECUTION_EXITED_ZERO
    assert report.entries[0].output_status == OUTPUT_NO_OUTPUT
    assert report.entries[0].status == COMPLETED_INVALID
    evidence = json.loads(
        (stage.stage_root / report.entries[0].atom_path / "rsutil_execution.json").read_text()
    )
    assert evidence["identity"]["exchange_stage_id"] == stage.stage_id
    assert evidence["identity"]["parent_entry_id"] == entry.entry_id


def test_collection_uses_canonical_exchange_parser_and_preserves_pair_provenance(tmp_path: Path) -> None:
    entry = _write_parent(tmp_path, entry_id="a" * 64, name="valid")
    _write_campaign(tmp_path, (entry,))
    stage = prepare_exchange_stage(tmp_path, tmp_path / "exchange-stage")
    job = stage.jobs[0]
    atom = stage.stage_root / job["atom_path"]
    manifest = json.loads(
        (stage.stage_root / job["exchange_manifest_path"]).read_text()
    )
    rows = []
    for pair in manifest["spec"]["pairs"]:
        rows.append(
            " ".join(
                str(value)
                for value in (
                    pair["rs_type_i"],
                    pair["rs_type_j"],
                    *pair["cartesian_displacement"],
                    -1.25,
                    pair["distance"],
                )
            )
        )
    (atom / "jijso.out").write_text("\n".join(rows) + "\n", encoding="utf-8")
    (atom / "rsutil_execution.json").write_text(
        json.dumps({
            "schema_version": "rsutil-execution-v1",
            "identity": {
                "exchange_stage_id": stage.stage_id,
                "campaign_id": "campaign-id",
                "parent_entry_id": job["parent_entry_id"],
                "exchange_manifest_hash": job["exchange_manifest_hash"],
                "job_identity": job["job_identity"],
                "array_index": "0",
            },
            "status": "success",
            "return_code": 0,
        }),
        encoding="utf-8",
    )

    report = status_exchange_stage(stage.stage_root)
    assert report.entries[0].status == "COMPLETED_VALID"
    collected = collect_exchange_stage(stage.stage_root)
    payload = json.loads(collected.aggregate_path.read_text())
    record = payload["entries"][0]["records"][0]
    assert record["key"]
    assert record["pair"]["rs_type_i"] is not None
    assert record["pair"]["species_i"] == "Fe"
