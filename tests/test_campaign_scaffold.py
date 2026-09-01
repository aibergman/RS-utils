"""Coverage for RS-LMTO source-bundle campaign scaffolding."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from rsutil.campaign import CampaignScaffoldError, load_campaign, scaffold_campaign, validate_campaign
from rsutil.core.hashing import hash_file

from .conftest import copy_example_files, run_rsutil


def _make_bundle(
    root: Path,
    relative: str,
    *,
    source_cif: bool = False,
    cif_files: tuple[str, ...] = (),
) -> Path:
    directory = root / relative
    directory.mkdir(parents=True)
    copy_example_files(directory, "bccFe", ["input.nml", "lattice.nml"])
    if source_cif:
        (directory / "source.cif").write_text("data-example\n", encoding="utf-8")
    for filename in cif_files:
        (directory / filename).write_text(f"data-{filename}\n", encoding="utf-8")
    return directory


def test_scaffold_writes_current_campaign_schema_with_relative_sources(tmp_path: Path) -> None:
    source_root = tmp_path / "sources"
    _make_bundle(source_root, "mp-13--Fe", source_cif=True)
    _make_bundle(source_root, "nested/mp-1007824--MnSi")
    output = tmp_path / "campaigns" / "campaign.json"

    result = scaffold_campaign(source_root, name="initial-materials", output=output)
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert payload["schema_version"] == "rsutil-campaign-v1"
    assert payload["name"] == "initial-materials"
    assert result.discovered == 2
    assert result.with_source_cif == 1
    assert all(not Path(path).is_absolute() for path in payload["structures"].values())
    assert all(entry["calculation"] == {"mode": "scf"} for entry in payload["entries"])
    assert all(entry["structure"] in payload["structures"] for entry in payload["entries"])
    fe_entry = next(entry for entry in payload["entries"] if entry["name"] == "mp-13--Fe")
    assert fe_entry["metadata"]["structure_source"] == {
        "kind": "cif",
        "file": "source.cif",
        "sha256": hash_file(source_root / "mp-13--Fe" / "source.cif"),
    }

    loaded = load_campaign(output)
    assert loaded.name == "initial-materials"
    assert len(loaded.entries) == 2
    assert loaded.entries[0].metadata == payload["entries"][0]["metadata"]


def test_scaffold_cif_hash_changes_campaign_entry_identity(tmp_path: Path) -> None:
    source_root = tmp_path / "sources"
    bundle = _make_bundle(source_root, "Fe", source_cif=True)
    output = tmp_path / "campaign.json"

    first = scaffold_campaign(source_root, name="materials", output=output)
    first_metadata = first.campaign.entries[0].metadata
    first_id = next(iter(first.validation.entry_ids.values()))

    (bundle / "source.cif").write_text("data-changed\n", encoding="utf-8")
    second = scaffold_campaign(source_root, name="materials", output=output)
    second_metadata = second.campaign.entries[0].metadata
    second_id = next(iter(second.validation.entry_ids.values()))

    assert first_metadata["structure_source"]["sha256"] != second_metadata["structure_source"]["sha256"]
    assert first_id != second_id
    assert second_metadata["structure_source"]["sha256"] == hash_file(bundle / "source.cif")


def test_scaffold_without_cif_has_no_source_metadata_and_relocates(tmp_path: Path) -> None:
    source_root = tmp_path / "sources"
    _make_bundle(source_root, "Fe")
    output = tmp_path / "campaign.json"

    scaffold_campaign(source_root, name="materials", output=output)
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["entries"][0]["metadata"] == {}
    assert payload["structures"]["Fe"] == "sources/Fe"

    relocated = tmp_path.parent / f"{tmp_path.name}-relocated"
    shutil.copytree(tmp_path, relocated)
    relocated_campaign = load_campaign(relocated / "campaign.json")
    assert validate_campaign(relocated_campaign, base_dir=relocated).valid


def test_scaffold_does_not_guess_among_ambiguous_cifs(tmp_path: Path) -> None:
    source_root = tmp_path / "sources"
    _make_bundle(source_root, "Fe", cif_files=("a.cif", "b.cif"))

    result = scaffold_campaign(source_root, name="materials", output=tmp_path / "campaign.json")

    assert result.ambiguous_cif == (("Fe", ("a.cif", "b.cif")),)
    assert result.campaign.entries[0].metadata == {}


def test_scaffold_cli_reports_summary(tmp_path: Path) -> None:
    source_root = tmp_path / "sources"
    _make_bundle(source_root, "Fe", source_cif=True)
    output = tmp_path / "campaign.json"

    completed = run_rsutil(
        tmp_path,
        "campaign",
        "scaffold",
        "sources",
        "--name",
        "initial-materials",
        "--output",
        "campaign.json",
    )

    assert completed.returncode == 0, completed.stderr
    assert "Discovered source bundles:  1" in completed.stdout
    assert "With source.cif:             1" in completed.stdout
    assert "Validation: PASS" in completed.stdout
    assert output.is_file()


def test_scaffold_fails_closed_on_incomplete_candidate(tmp_path: Path) -> None:
    source_root = tmp_path / "sources"
    _make_bundle(source_root, "good")
    incomplete = source_root / "broken-Fe"
    incomplete.mkdir(parents=True)
    (incomplete / "input.nml").write_text("&atoms\n/\n", encoding="utf-8")
    output = tmp_path / "campaign.json"

    with pytest.raises(CampaignScaffoldError, match=r"broken-Fe/ missing lattice\.nml"):
        scaffold_campaign(source_root, name="initial-materials", output=output)
    assert not output.exists()


def test_scaffold_rejects_sanitized_alias_collision(tmp_path: Path) -> None:
    source_root = tmp_path / "sources"
    _make_bundle(source_root, "A B")
    _make_bundle(source_root, "A-B")

    with pytest.raises(CampaignScaffoldError, match="source alias collision"):
        scaffold_campaign(source_root, name="collision", output=tmp_path / "campaign.json")
