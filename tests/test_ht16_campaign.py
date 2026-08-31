"""HT-16 deterministic campaign specification and generation coverage."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from rsutil.campaign import Campaign, CampaignCollisionError, CampaignGenerator
from rsutil.core import OperationStatus, SetupResult, ValidationResult
from rsutil.setup import SetupService

from .conftest import copy_example_files


class _FakeRunner:
    def run(self, lattice_path, input_path, spec, *, workspace):
        workspace = Path(workspace)
        workspace.mkdir(parents=True, exist_ok=True)
        outputs = {
            "lattice_new.nml": "&lattice\n ntot = 1\n nbas = 1\n nrec = 1\n a = 1,0,0, 0,1,0, 0,0,1\n crd = 0,0,0\n izp = 26\n no = 1\n iu = 1\n ib = 1\n irec = 1\n nbulk_bulk = 1\n/\n",
            "input_new.nml": "&atoms\n label(1) = 'Fe'\n/\n",
            "radii.out": "1\n1.0\n",
        }
        for name, content in outputs.items():
            (workspace / name).write_text(content, encoding="utf-8")
        return SetupResult(
            status=OperationStatus.SUCCESS,
            candidate_files=tuple(outputs),
            semantic_validation=ValidationResult(valid=True),
        )


def _source(root: Path, example: str, destination: str) -> None:
    target = root / destination
    target.mkdir()
    names = ["Dist/lattice.nml", "Dist/input.nml"] if example == "bccFe" else ["lattice.nml", "input.nml"]
    copy_example_files(target, example, names)


def _campaign(source_root: Path, *, second_setup: dict | None = None) -> Campaign:
    entries = [
        {"name": "reference", "source": "bcc", "setup": {"seed": 41}, "calculation": {"mode": "scf"}},
        {"name": "second", "source": "fen", "setup": second_setup or {"seed": 42}, "calculation": {"mode": "scf"}},
    ]
    return Campaign.from_dict({
        "name": "pilot",
        "structures": {"bcc": "bcc", "fen": "fen"},
        "entries": entries,
    })


def test_campaign_generation_is_idempotent_and_relocatable(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    _source(tmp_path, "FeN", "fen")
    campaign = _campaign(tmp_path)
    generator = CampaignGenerator(setup_service=SetupService(runner=_FakeRunner()))
    destination = tmp_path / "campaign"

    first = generator.generate(campaign, destination, base_dir=tmp_path)
    first_bytes = {path: path.read_bytes() for path in destination.rglob("*") if path.is_file()}
    second = generator.generate(campaign, destination, base_dir=tmp_path)

    assert [entry.entry_id for entry in first.entries] == [entry.entry_id for entry in second.entries]
    assert first_bytes == {path: path.read_bytes() for path in destination.rglob("*") if path.is_file()}
    assert len(list(destination.glob("runs/*"))) == 2
    assert all(not Path(entry.run_path).is_absolute() for entry in second.entries)

    relocated = tmp_path / "relocated"
    shutil.copytree(destination, relocated)
    index = json.loads((relocated / "campaign_index.json").read_text(encoding="utf-8"))
    assert all(not str(value).startswith("/") for entry in index["entries"] for value in (entry["run_path"], entry["manifest_path"]))
    assert all((relocated / entry["manifest_path"]).is_file() for entry in index["entries"])


def test_only_changed_scientific_entry_identity_changes(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    _source(tmp_path, "FeN", "fen")
    original = _campaign(tmp_path)
    changed = _campaign(tmp_path, second_setup={"seed": 42, "lmx": 3})
    original_validation = CampaignGenerator().validate(original, base_dir=tmp_path)
    changed_validation = CampaignGenerator().validate(changed, base_dir=tmp_path)

    assert original_validation.entry_ids["reference"] == changed_validation.entry_ids["reference"]
    assert original_validation.entry_ids["second"] != changed_validation.entry_ids["second"]


def test_invalid_entry_fails_before_destination_mutation(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    campaign = Campaign.from_dict({
        "name": "invalid",
        "entries": [
            {"name": "good", "source": "bcc"},
            {"name": "bad", "source": "does-not-exist"},
        ],
    })
    destination = tmp_path / "not-created"

    with pytest.raises(Exception, match="preflight validation failed"):
        CampaignGenerator().generate(campaign, destination, base_dir=tmp_path)
    assert not destination.exists()


def test_modified_spec_under_explicit_same_name_path_is_a_collision(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    first = Campaign.from_dict({
        "name": "collision",
        "destination_naming": "runs/shared",
        "entries": [{"name": "case", "source": "bcc", "setup": {"seed": 1}}],
    })
    second = Campaign.from_dict({
        "name": "collision",
        "destination_naming": "runs/shared",
        "entries": [{"name": "case", "source": "bcc", "setup": {"seed": 2}}],
    })
    generator = CampaignGenerator(setup_service=SetupService(runner=_FakeRunner()))
    generator.generate(first, tmp_path / "output", base_dir=tmp_path)
    with pytest.raises(CampaignCollisionError):
        generator.generate(second, tmp_path / "output", base_dir=tmp_path)
