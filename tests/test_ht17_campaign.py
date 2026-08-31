"""HT-17 Slurm rendering, filesystem status, and collection coverage."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

from rsutil.campaign import (
    CAMPAIGN_RESULTS_SCHEMA_VERSION,
    COMPLETED_INVALID,
    COMPLETED_VALID,
    FAILED,
    GENERATED,
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
        "valid case": COMPLETED_VALID,
        "malformed": COMPLETED_INVALID,
        "generated": GENERATED,
        "failed": FAILED,
    }
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
    assert valid["result_sidecar"] == f"collected/results/{index.entries[0].entry_id}.json"
    assert (tmp_path / valid["result_sidecar"]).is_file()
    assert malformed["parse_error"]["code"] == "malformed_result"
    assert result.csv_path is not None and result.csv_path.is_file()
    assert result.jsonl_path is not None and len(result.jsonl_path.read_text(encoding="utf-8").splitlines()) == 4
