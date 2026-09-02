"""HT-OPS-01 resource planning and heterogeneous Slurm coverage."""

from __future__ import annotations

import json
import os
import shlex
import subprocess
from pathlib import Path

from rsutil.campaign import CampaignEntry, CampaignIndex, CampaignIndexEntry, render_slurm_array


def _index(root: Path) -> CampaignIndex:
    entries = tuple(
        CampaignIndexEntry(
            entry_id=letter * 64,
            name=name,
            run_path=f"runs/{name}",
            manifest_path=f"runs/{name}/run_manifest.json",
            manifest_hash=f"manifest-{letter}",
            resources=resources,
        )
        for letter, name, resources in (
            ("a", "large path", {"ntasks": 8, "time": "02:00:00"}),
            ("b", "small", {}),
            ("c", "medium", {"ntasks": 4, "time": "01:00:00"}),
            ("d", "medium copy", {"ntasks": 4, "time": "01:00:00"}),
            ("e", "small copy", {"ntasks": 2, "time": "00:30:00"}),
        )
    )
    index = CampaignIndex("campaign-id", "resource campaign", "spec-hash", entries)
    for entry in entries:
        run = root / entry.run_path
        run.mkdir(parents=True)
        (run / "run_manifest.json").write_text(
            json.dumps({"campaign": {"entry_id": entry.entry_id, "manifest_hash": entry.manifest_hash}}),
            encoding="utf-8",
        )
    (root / "campaign_index.json").write_text(json.dumps(index.to_dict()), encoding="utf-8")
    return index


def test_resource_settings_are_not_scientific_identity(tmp_path: Path) -> None:
    first = CampaignEntry("case", source="source", resources={"ntasks": 2})
    second = CampaignEntry("case", source="source", resources={"ntasks": 8, "time": "02:00:00"})

    assert first.identity_hash(source_hashes={"input.nml": "same"}) == second.identity_hash(
        source_hashes={"input.nml": "same"}
    )


def test_heterogeneous_rendering_groups_deterministically_and_inherits_defaults(tmp_path: Path) -> None:
    index = _index(tmp_path)
    fake = tmp_path / "fake executable.sh"
    fake.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    fake.chmod(0o755)

    first = render_slurm_array(
        index,
        campaign_root=tmp_path,
        output_dir=tmp_path / "slurm",
        executable=fake,
        time_limit="00:10:00",
        resources={"cpus_per_task": 1, "memory": "4G"},
    )
    second = render_slurm_array(
        index,
        campaign_root=tmp_path,
        output_dir=tmp_path / "slurm-second",
        executable=fake,
        time_limit="00:10:00",
        resources={"cpus_per_task": 1, "mem": "4G"},
    )

    assert first.resource_aware is True
    assert first.class_count == 4
    assert [path.name for path in sorted((tmp_path / "slurm").glob("class-*/"))] == [
        "class-001", "class-002", "class-003", "class-004"
    ]
    plan = json.loads((tmp_path / "slurm" / "resource_plan.json").read_text(encoding="utf-8"))
    assert [item["resource_class"] for item in plan["classes"]] == [
        "class-001", "class-002", "class-003", "class-004"
    ]
    by_name = {item["entry_id"]: item for item in plan["entries"]}
    assert by_name["b" * 64]["ntasks"] is None
    assert by_name["b" * 64]["cpus_per_task"] == 1
    assert by_name["b" * 64]["memory"] == "4G"
    assert by_name["b" * 64]["time"] == "00:10:00"
    assert by_name["a" * 64]["ntasks"] == 8
    assert by_name["a" * 64]["time"] == "02:00:00"
    assert by_name["c" * 64]["resource_class"] == by_name["d" * 64]["resource_class"]

    submit_all = (tmp_path / "slurm" / "submit_all.sh").read_text(encoding="utf-8")
    submit_lines = [line for line in submit_all.splitlines() if line.startswith("sbatch ")]
    assert len(submit_lines) == first.class_count
    assert {Path(shlex.split(line)[1]).parent.name for line in submit_lines} == {
        f"class-{number:03d}" for number in range(1, first.class_count + 1)
    }
    assert json.loads((tmp_path / "slurm-second" / "resource_plan.json").read_text(encoding="utf-8")) == plan


def test_resource_aware_homogeneous_campaign_keeps_one_class_and_quotes_paths(tmp_path: Path) -> None:
    index = _index(tmp_path)
    record = tmp_path / "selected result.txt"
    fake = tmp_path / "fake executable.sh"
    fake.write_text(
        "#!/usr/bin/env bash\n"
        f"printf '%s|%s|%s\\n' \"$RSUTIL_CAMPAIGN_ARRAY_INDEX\" \"$RSUTIL_CAMPAIGN_ENTRY_ID\" \"$PWD\" > {shlex.quote(str(record))}\n",
        encoding="utf-8",
    )
    fake.chmod(0o755)

    rendered = render_slurm_array(
        index,
        campaign_root=tmp_path,
        output_dir=tmp_path / "slurm",
        executable=fake,
        resource_aware=True,
        resources={"ntasks": 2},
    )
    assert rendered.class_count == 4
    assert "#SBATCH --ntasks=2" in (tmp_path / "slurm" / "class-002" / "submit_array.sbatch").read_text()
    environment = os.environ.copy()
    environment["SLURM_ARRAY_TASK_ID"] = "0"
    completed = subprocess.run(["bash", str(rendered.script_path)], env=environment, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert record.read_text(encoding="utf-8").startswith("0|" + "b" * 64 + "|")

