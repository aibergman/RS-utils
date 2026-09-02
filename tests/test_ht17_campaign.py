"""HT-17 Slurm rendering, filesystem status, and collection coverage."""

from __future__ import annotations

import csv
import json
import os
import subprocess
from pathlib import Path

from rsutil.campaign import (
    CAMPAIGN_RESULTS_SCHEMA_VERSION,
    COMPLETED_INVALID,
    COMPLETED_NOT_CONVERGED,
    COMPLETED_UNKNOWN,
    COMPLETED_VALID,
    EXECUTION_EXITED_NONZERO,
    EXECUTION_EXITED_ZERO,
    EXECUTION_RUNNING_OR_UNKNOWN,
    EXECUTION_SIGNALED,
    FAILED,
    GENERATED,
    OUTPUT_NO_OUTPUT,
    OUTPUT_PARSE_INVALID,
    OUTPUT_PARSE_VALID,
    SCIENTIFIC_CONVERGED,
    SCIENTIFIC_NOT_CONVERGED,
    SCIENTIFIC_UNKNOWN,
    CampaignIndex,
    CampaignIndexEntry,
    collect_campaign,
    render_slurm_array,
    status_campaign,
)


def _make_index(root: Path) -> CampaignIndex:
    entries = (
        CampaignIndexEntry("a" * 64, "valid case", "runs/valid case", "runs/valid case/run_manifest.json", "manifest-a"),
        CampaignIndexEntry("b" * 64, "malformed", "runs/malformed;case", "runs/malformed;case/run_manifest.json", "manifest-b"),
        CampaignIndexEntry("c" * 64, "generated", "runs/generated", "runs/generated/run_manifest.json", "manifest-c"),
        CampaignIndexEntry("d" * 64, "failed", "runs/failed", "runs/failed/run_manifest.json", "manifest-d"),
    )
    index = CampaignIndex("campaign-id", "spaces & quotes", "spec-hash", entries)
    for entry in entries:
        run_dir = root / entry.run_path
        run_dir.mkdir(parents=True)
        (run_dir / "run_manifest.json").write_text(
            json.dumps({
                "schema_version": "rsutil-run-v1",
                "campaign": {
                    "campaign_name": index.campaign_name,
                    "entry_name": entry.name,
                    "entry_id": entry.entry_id,
                    "manifest_hash": entry.manifest_hash,
                },
            }),
            encoding="utf-8",
        )
    valid = root / entries[0].run_path
    (valid / "report.out").write_text("Fermi energy: -0.125\n", encoding="utf-8")
    (valid / "out.scf").write_text("Converged!\n", encoding="utf-8")
    malformed = root / entries[1].run_path
    (malformed / "report.out").write_text("Fermi energy: not-a-number\n", encoding="utf-8")
    (root / entries[3].run_path / "FAILED").write_text("fake executable exit 1\n", encoding="utf-8")
    (root / "campaign_index.json").write_text(json.dumps(index.to_dict()), encoding="utf-8")
    return index


def test_render_slurm_is_explicit_and_shell_safe(tmp_path: Path) -> None:
    index = _make_index(tmp_path)
    fake = tmp_path / "fake rs executable.sh"
    record = tmp_path / "selected entry.txt"
    fake.write_text(
        "#!/usr/bin/env bash\n"
        f"printf '%s|%s|%s\\n' \"$RSUTIL_CAMPAIGN_ARRAY_INDEX\" \"$RSUTIL_CAMPAIGN_ENTRY_ID\" \"$PWD\" > {json.dumps(str(record))}\n",
        encoding="utf-8",
    )
    fake.chmod(fake.stat().st_mode | 0o111)

    rendered = render_slurm_array(index, campaign_root=tmp_path, executable=str(fake))
    run_list = rendered.run_list_path.read_text(encoding="utf-8")
    script = rendered.script_path.read_text(encoding="utf-8")

    assert run_list.splitlines()[0] == "array_index\tentry_id\trun_path\tmanifest_path"
    assert "1\t" + "b" * 64 + "\truns/malformed;case" in run_list
    assert "#SBATCH --array=0-3" in script
    assert "fake rs executable.sh" in script
    assert "module load" not in script

    environment = os.environ.copy()
    environment["SLURM_ARRAY_TASK_ID"] = "1"
    completed = subprocess.run(["bash", str(rendered.script_path)], env=environment, capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr
    assert record.read_text(encoding="utf-8").startswith("1|" + "b" * 64 + "|")
    assert record.read_text(encoding="utf-8").rstrip().endswith("runs/malformed;case")


def test_status_uses_filesystem_evidence_without_scheduler_claims(tmp_path: Path) -> None:
    _make_index(tmp_path)
    report = status_campaign(tmp_path)
    statuses = {entry.name: entry.status for entry in report.entries}
    assert statuses == {
        "valid case": "UNKNOWN",
        "malformed": COMPLETED_INVALID,
        "generated": GENERATED,
        "failed": FAILED,
    }
    valid = next(entry for entry in report.entries if entry.name == "valid case")
    assert valid.execution_status == EXECUTION_RUNNING_OR_UNKNOWN
    assert valid.output_status == OUTPUT_PARSE_VALID
    assert valid.scientific_status == SCIENTIFIC_CONVERGED
    assert report.scheduler_state == "not_queried"


def test_collection_continues_and_preserves_hashes_and_sidecars(tmp_path: Path) -> None:
    index = _make_index(tmp_path)
    result = collect_campaign(tmp_path, csv_output=tmp_path / "summary.csv", jsonl_output=tmp_path / "summary.jsonl")
    aggregate = json.loads(result.aggregate_path.read_text(encoding="utf-8"))
    assert aggregate["schema_version"] == CAMPAIGN_RESULTS_SCHEMA_VERSION
    assert len(aggregate["entries"]) == 4
    valid = next(entry for entry in aggregate["entries"] if entry["entry_id"] == index.entries[0].entry_id)
    malformed = next(entry for entry in aggregate["entries"] if entry["entry_id"] == index.entries[1].entry_id)
    assert valid["manifest_hash"] == index.entries[0].manifest_hash
    assert valid["manifest_content_hash"]
    assert valid["execution_status"] == EXECUTION_RUNNING_OR_UNKNOWN
    assert valid["output_status"] == OUTPUT_PARSE_VALID
    assert valid["scientific_status"] == SCIENTIFIC_CONVERGED
    assert valid["result_sidecar"] == f"collected/results/{index.entries[0].entry_id}.json"
    assert (tmp_path / valid["result_sidecar"]).is_file()
    assert malformed["parse_error"]["code"] == "malformed_result"
    assert result.csv_path is not None and result.csv_path.is_file()
    assert result.jsonl_path is not None and len(result.jsonl_path.read_text(encoding="utf-8").splitlines()) == 4


def test_collection_exports_vector_magnetic_scalars_without_classifying_entries(tmp_path: Path) -> None:
    index = _make_index(tmp_path)
    valid = tmp_path / index.entries[0].run_path
    (valid / "report.out").write_text(
        "Total spin moment: 1.0 2.0 2.0\n"
        "Total orbital moment: 0.0 0.0 0.5\n"
        "Fermi energy: -0.125\n",
        encoding="utf-8",
    )
    collection = collect_campaign(
        tmp_path,
        csv_output=tmp_path / "summary.csv",
        jsonl_output=tmp_path / "summary.jsonl",
    )

    entry = collection.entries[0]
    summary = entry["summary"]
    assert summary["total_spin_moment_x"] == 1.0
    assert summary["total_spin_moment_y"] == 2.0
    assert summary["total_spin_moment_z"] == 2.0
    assert summary["total_spin_moment_norm"] == 3.0
    assert summary["total_orbital_moment_norm"] == 0.5
    assert "magnetic" not in summary
    assert "nonmagnetic" not in summary

    jsonl_entry = json.loads(collection.jsonl_path.read_text(encoding="utf-8").splitlines()[0])
    assert jsonl_entry["summary"]["total_spin_moment_norm"] == 3.0
    with collection.csv_path.open(newline="", encoding="utf-8") as stream:
        row = next(csv.DictReader(stream))
    assert float(row["total_spin_moment_norm"]) == 3.0
    assert "magnetic" not in row
    assert "nonmagnetic" not in row


def test_collection_keeps_missing_magnetic_scalars_empty(tmp_path: Path) -> None:
    index = _make_index(tmp_path)
    collection = collect_campaign(tmp_path, csv_output=tmp_path / "summary.csv")
    summary = collection.entries[0]["summary"]
    assert summary["total_spin_moment_x"] is None
    assert summary["total_spin_moment_y"] is None
    assert summary["total_spin_moment_z"] is None
    assert summary["total_spin_moment_norm"] is None
    with collection.csv_path.open(newline="", encoding="utf-8") as stream:
        row = next(csv.DictReader(stream))
    assert row["total_spin_moment_norm"] == ""


def _single_index(root: Path, name: str = "case") -> tuple[CampaignIndex, Path]:
    entry = CampaignIndexEntry("e" * 64, name, f"runs/{name}", f"runs/{name}/run_manifest.json", "manifest-e")
    index = CampaignIndex("campaign-id", "campaign", "spec-hash", (entry,))
    run = root / entry.run_path
    run.mkdir(parents=True)
    (run / "run_manifest.json").write_text(json.dumps({
        "schema_version": "rsutil-run-v1",
        "campaign": {"entry_id": entry.entry_id, "manifest_hash": entry.manifest_hash},
    }), encoding="utf-8")
    (root / "campaign_index.json").write_text(json.dumps(index.to_dict()), encoding="utf-8")
    return index, run


def _execute_fake(
    root: Path,
    *,
    report: str | None,
    log: str | None,
    exit_code: int,
) -> tuple[CampaignIndex, Path]:
    index, run = _single_index(root)
    fake = root / "fake.sh"
    lines = ["#!/usr/bin/env bash"]
    if report is not None:
        lines.append(f"printf '%b' {json.dumps(report)} > report.out")
    if log is not None:
        lines.append(f"printf '%b' {json.dumps(log)} > out.scf")
    lines.append(f"exit {exit_code}")
    fake.write_text("\n".join(lines) + "\n", encoding="utf-8")
    fake.chmod(0o755)
    rendered = render_slurm_array(index, campaign_root=root, executable=fake)
    environment = os.environ.copy()
    environment["SLURM_ARRAY_TASK_ID"] = "0"
    completed = subprocess.run(["bash", str(rendered.script_path)], env=environment, capture_output=True, text=True, check=False)
    assert completed.returncode == exit_code, completed.stderr
    return index, run


def test_generated_runner_records_success_and_converged_axes(tmp_path: Path) -> None:
    index, run = _execute_fake(tmp_path, report="Fermi energy: -0.125\n", log="Converged!\n", exit_code=0)
    evidence = json.loads((run / "rsutil_execution.json").read_text(encoding="utf-8"))
    assert evidence["schema_version"] == "rsutil-execution-v1"
    assert evidence["status"] == "success"
    assert evidence["return_code"] == 0
    assert evidence["identity"]["entry_id"] == index.entries[0].entry_id
    entry = status_campaign(tmp_path).entries[0]
    assert entry.status == COMPLETED_VALID
    assert entry.execution_status == EXECUTION_EXITED_ZERO
    assert entry.output_status == OUTPUT_PARSE_VALID
    assert entry.scientific_status == SCIENTIFIC_CONVERGED


def test_nonconverged_report_is_not_a_scientific_success(tmp_path: Path) -> None:
    _execute_fake(tmp_path, report="Fermi energy: -0.125\nNot converged!\n", log="Not converged!\n", exit_code=0)
    entry = status_campaign(tmp_path).entries[0]
    assert entry.status == COMPLETED_NOT_CONVERGED
    assert entry.execution_status == EXECUTION_EXITED_ZERO
    assert entry.output_status == OUTPUT_PARSE_VALID
    assert entry.scientific_status == SCIENTIFIC_NOT_CONVERGED
    collection = collect_campaign(tmp_path)
    assert not collection.successful
    assert len(collection.failures) == 1
    assert collection.entries[0]["scientific_status"] == SCIENTIFIC_NOT_CONVERGED


def test_zero_exit_with_unknown_convergence_remains_explicit(tmp_path: Path) -> None:
    _execute_fake(tmp_path, report="Fermi energy: -0.125\n", log=None, exit_code=0)
    entry = status_campaign(tmp_path).entries[0]
    assert entry.status == COMPLETED_UNKNOWN
    assert entry.execution_status == EXECUTION_EXITED_ZERO
    assert entry.output_status == OUTPUT_PARSE_VALID
    assert entry.scientific_status == SCIENTIFIC_UNKNOWN


def test_zero_exit_with_malformed_report_is_parser_failure(tmp_path: Path) -> None:
    _execute_fake(tmp_path, report="Fermi energy: not-a-number\n", log=None, exit_code=0)
    entry = status_campaign(tmp_path).entries[0]
    assert entry.status == COMPLETED_INVALID
    assert entry.execution_status == EXECUTION_EXITED_ZERO
    assert entry.output_status == OUTPUT_PARSE_INVALID
    assert entry.scientific_status == SCIENTIFIC_UNKNOWN


def test_nonzero_exit_without_report_is_persisted_as_execution_failure(tmp_path: Path) -> None:
    _execute_fake(tmp_path, report=None, log=None, exit_code=9)
    entry = status_campaign(tmp_path).entries[0]
    assert entry.status == FAILED
    assert entry.execution_status == EXECUTION_EXITED_NONZERO
    assert entry.output_status == OUTPUT_NO_OUTPUT


def test_signaled_execution_is_distinguished_from_nonzero_exit(tmp_path: Path) -> None:
    index, _run = _single_index(tmp_path)
    fake = tmp_path / "signaled.sh"
    fake.write_text("#!/usr/bin/env bash\nkill -TERM $$\n", encoding="utf-8")
    fake.chmod(0o755)
    rendered = render_slurm_array(index, campaign_root=tmp_path, executable=fake)
    environment = os.environ.copy(); environment["SLURM_ARRAY_TASK_ID"] = "0"
    completed = subprocess.run(["bash", str(rendered.script_path)], env=environment, capture_output=True, check=False)
    assert completed.returncode != 0
    entry = status_campaign(tmp_path).entries[0]
    assert entry.status == FAILED
    assert entry.execution_status == EXECUTION_SIGNALED
    assert entry.execution_evidence["signal"] == 15


def test_nonzero_exit_cannot_be_hidden_by_stale_parseable_report(tmp_path: Path) -> None:
    index, run = _single_index(tmp_path)
    (run / "report.out").write_text("Fermi energy: -0.125\n", encoding="utf-8")
    (run / "out.scf").write_text("Converged!\n", encoding="utf-8")
    fake = tmp_path / "fails.sh"
    fake.write_text("#!/usr/bin/env bash\nexit 4\n", encoding="utf-8")
    fake.chmod(0o755)
    rendered = render_slurm_array(index, campaign_root=tmp_path, executable=fake)
    environment = os.environ.copy(); environment["SLURM_ARRAY_TASK_ID"] = "0"
    assert subprocess.run(["bash", str(rendered.script_path)], env=environment, capture_output=True, check=False).returncode == 4
    entry = status_campaign(tmp_path).entries[0]
    assert entry.status == FAILED
    assert entry.execution_status == EXECUTION_EXITED_NONZERO
    assert entry.output_status == OUTPUT_PARSE_VALID
    assert entry.scientific_status == SCIENTIFIC_CONVERGED


def test_legacy_report_without_execution_evidence_is_unknown_not_success(tmp_path: Path) -> None:
    _index, run = _single_index(tmp_path)
    (run / "report.out").write_text("Fermi energy: -0.125\n", encoding="utf-8")
    (run / "out.scf").write_text("Converged!\n", encoding="utf-8")
    entry = status_campaign(tmp_path).entries[0]
    assert entry.status == "UNKNOWN"
    assert entry.execution_status == EXECUTION_RUNNING_OR_UNKNOWN
    assert entry.output_status == OUTPUT_PARSE_VALID
    assert entry.scientific_status == SCIENTIFIC_CONVERGED


def test_wrong_execution_identity_is_not_accepted(tmp_path: Path) -> None:
    _index, run = _single_index(tmp_path)
    (run / "report.out").write_text("Fermi energy: -0.125\n", encoding="utf-8")
    (run / "out.scf").write_text("Converged!\n", encoding="utf-8")
    (run / "rsutil_execution.json").write_text(json.dumps({
        "schema_version": "rsutil-execution-v1",
        "identity": {"entry_id": "wrong", "manifest_hash": "manifest-e"},
        "status": "success", "return_code": 0,
    }), encoding="utf-8")
    entry = status_campaign(tmp_path).entries[0]
    assert entry.status == "UNKNOWN"
    assert entry.execution_status == EXECUTION_RUNNING_OR_UNKNOWN
    assert entry.parser_error["code"] == "execution_identity_mismatch"
