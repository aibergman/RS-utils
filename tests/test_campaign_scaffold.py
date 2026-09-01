"""Coverage for RS-LMTO source-bundle campaign scaffolding."""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path

import pytest

from rsutil.campaign import CampaignGenerator, CampaignScaffoldError, load_campaign, scaffold_campaign, validate_campaign
from rsutil.core import OperationStatus, SetupResult, ValidationResult
from rsutil.core.hashing import hash_file
from rsutil.setup import SetupService

from .conftest import copy_example_files, run_rsutil


class _FakeRunner:
    """Keep the lifecycle test focused on campaign plumbing, not Fortran."""

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
    assert all(entry["metadata"] == {} for entry in payload["entries"])
    assert payload["metadata"]["structure_sources"] == {
        "mp-13--Fe": {
            "kind": "cif",
            "file": "source.cif",
            "sha256": hash_file(source_root / "mp-13--Fe" / "source.cif"),
        }
    }

    loaded = load_campaign(output)
    assert loaded.name == "initial-materials"
    assert len(loaded.entries) == 2
    assert loaded.metadata == payload["metadata"]
    assert loaded.entries[0].metadata == {}


def test_scaffold_cif_snapshot_changes_campaign_provenance_not_entry_identity(tmp_path: Path) -> None:
    source_root = tmp_path / "sources"
    bundle = _make_bundle(source_root, "Fe", source_cif=True)
    first_output = tmp_path / "campaign-first.json"
    second_output = tmp_path / "campaign-second.json"

    first = scaffold_campaign(source_root, name="materials", output=first_output)
    first_metadata = first.campaign.metadata
    first_id = next(iter(first.validation.entry_ids.values()))

    (bundle / "source.cif").write_text("data-changed\n", encoding="utf-8")
    second = scaffold_campaign(source_root, name="materials", output=second_output)
    second_metadata = second.campaign.metadata
    second_id = next(iter(second.validation.entry_ids.values()))

    assert first_metadata["structure_sources"]["Fe"]["sha256"] != second_metadata["structure_sources"]["Fe"]["sha256"]
    assert first_id == second_id
    assert second_metadata["structure_sources"]["Fe"]["sha256"] == hash_file(bundle / "source.cif")


def test_scaffold_without_cif_has_no_source_metadata_and_relocates(tmp_path: Path) -> None:
    source_root = tmp_path / "sources"
    _make_bundle(source_root, "Fe")
    output = tmp_path / "campaign.json"

    scaffold_campaign(source_root, name="materials", output=output)
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["entries"][0]["metadata"] == {}
    assert payload["metadata"] == {"structure_sources": {}}
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
    assert result.campaign.metadata == {"structure_sources": {}}


@pytest.mark.parametrize(
    ("filename", "old", "new"),
    [
        ("lattice.nml", "ct(1) = 3.3", "ct(1) = 3.4"),
        ("input.nml", "alat = 2.66000", "alat = 2.76000"),
    ],
)
def test_scaffold_scientific_source_changes_entry_identity(
    tmp_path: Path, filename: str, old: str, new: str
) -> None:
    source_root = tmp_path / "sources"
    bundle = _make_bundle(source_root, "Fe")
    first = scaffold_campaign(source_root, name="materials", output=tmp_path / "first.json")
    first_id = next(iter(first.validation.entry_ids.values()))

    path = bundle / filename
    path.write_text(path.read_text(encoding="utf-8").replace(old, new), encoding="utf-8")
    second = scaffold_campaign(source_root, name="materials", output=tmp_path / "second.json")
    second_id = next(iter(second.validation.entry_ids.values()))

    assert first_id != second_id


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


def _source_tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hash_file(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def test_scaffold_filters_are_deterministic_and_dry_run_is_non_mutating(tmp_path: Path) -> None:
    source_root = tmp_path / "sources"
    _make_bundle(source_root, "Fe", source_cif=True)
    _make_bundle(source_root, "test-Fe")
    _make_bundle(source_root, "MnSi")
    before = _source_tree_hashes(source_root)
    output = tmp_path / "campaign.json"

    first = scaffold_campaign(
        source_root,
        name="materials",
        output=output,
        include=("*Fe*",),
        exclude=("*test*",),
    )
    first_bytes = output.read_bytes()
    second = scaffold_campaign(
        source_root,
        name="materials",
        output=output,
        include=("*Fe*",),
        exclude=("*test*",),
    )

    assert first_bytes == output.read_bytes()
    assert first.discovered == second.discovered == 1
    assert [entry.name for entry in first.campaign.entries] == ["Fe"]
    assert _source_tree_hashes(source_root) == before

    dry_output = tmp_path / "dry-run.json"
    completed = run_rsutil(
        tmp_path,
        "campaign",
        "scaffold",
        "sources",
        "--name",
        "materials",
        "--output",
        "dry-run.json",
        "--include",
        "*Fe*",
        "--exclude",
        "*test*",
        "--dry-run",
    )
    assert completed.returncode == 0, completed.stderr
    assert "Dry run: no campaign JSON was written." in completed.stdout
    assert not dry_output.exists()
    assert _source_tree_hashes(source_root) == before


def test_scaffold_validate_create_lifecycle_preserves_sources_and_records_inputs(tmp_path: Path) -> None:
    source_root = tmp_path / "scaffold_sources"
    repository_sources = Path(__file__).resolve().parents[1] / "examples" / "campaigns" / "scaffold_sources"
    shutil.copytree(repository_sources, source_root)
    before = _source_tree_hashes(source_root)
    spec_path = tmp_path / "scaffolded.json"

    scaffold_campaign(source_root, name="scaffolded-pilot", output=spec_path)
    campaign = load_campaign(spec_path)
    validation = validate_campaign(campaign, base_dir=spec_path.parent)
    assert validation.valid
    assert set(validation.entry_ids) == {"FeCo", "FeN", "bcc_Fe"}

    destination = tmp_path / "pilot-runs"
    index = CampaignGenerator(
        setup_service=SetupService(runner=_FakeRunner())
    ).generate(campaign, destination, base_dir=spec_path.parent)
    assert len(index.entries) == 3
    assert {entry.name for entry in index.entries} == set(validation.entry_ids)
    for entry in index.entries:
        manifest = json.loads((destination / entry.manifest_path).read_text(encoding="utf-8"))
        source = manifest["campaign"]["source"]
        source_inputs = manifest["campaign"]["source_inputs"]
        assert source in {"FeCo", "FeN", "bcc_Fe"}
        assert source_inputs["input.nml"]["path"] == f"scaffold_sources/{source}/input.nml"
        assert source_inputs["lattice.nml"]["path"] == f"scaffold_sources/{source}/lattice.nml"
        assert source_inputs["input.nml"]["sha256"] == hash_file(source_root / source / "input.nml")
        assert source_inputs["lattice.nml"]["sha256"] == hash_file(source_root / source / "lattice.nml")
    assert _source_tree_hashes(source_root) == before


def test_scaffold_cif_is_optional_after_snapshot_for_validate_and_create(tmp_path: Path) -> None:
    source_root = tmp_path / "sources"
    bundle = _make_bundle(source_root, "Fe", source_cif=True)
    spec_path = tmp_path / "campaign.json"

    scaffold_campaign(source_root, name="materials", output=spec_path)
    bundle.joinpath("source.cif").unlink()
    campaign = load_campaign(spec_path)

    assert validate_campaign(campaign, base_dir=tmp_path).valid
    index = CampaignGenerator(setup_service=SetupService(runner=_FakeRunner())).generate(
        campaign, tmp_path / "runs", base_dir=tmp_path
    )
    assert len(index.entries) == 1


def test_scaffold_handles_large_source_collection_deterministically(tmp_path: Path) -> None:
    source_root = tmp_path / "many-sources"
    source_root.mkdir()
    fixture_root = Path(__file__).resolve().parents[1] / "examples" / "bccFe" / "Dist"
    for index in range(120):
        bundle = source_root / f"material-{index:03d}"
        bundle.mkdir()
        shutil.copy2(fixture_root / "input.nml", bundle / "input.nml")
        shutil.copy2(fixture_root / "lattice.nml", bundle / "lattice.nml")

    first_output = tmp_path / "many-first.json"
    second_output = tmp_path / "many-second.json"
    started = time.perf_counter()
    first = scaffold_campaign(source_root, name="many", output=first_output)
    elapsed = time.perf_counter() - started
    second = scaffold_campaign(source_root, name="many", output=second_output)

    assert first.discovered == first.complete == 120
    assert len(first.campaign.entries) == 120
    assert first_output.read_bytes() == second_output.read_bytes()
    assert elapsed < 30.0
