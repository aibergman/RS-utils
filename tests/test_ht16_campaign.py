"""HT-16 deterministic campaign specification and generation coverage."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from rsutil.campaign import Campaign, CampaignCollisionError, CampaignError, CampaignGenerator
from rsutil.core import OperationStatus, SetupResult, ValidationResult
from rsutil.core.hashing import hash_file
from rsutil.campaign.generation import _campaign_integrity_hash
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


def test_setup_preset_is_accepted_and_identity_bearing(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    base_data = {
        "name": "preset-identity",
        "structures": {"bcc": "bcc"},
        "entries": [{"name": "case", "structure": "bcc"}],
    }
    base = Campaign.from_dict(base_data)
    with_ks = Campaign.from_dict({**base_data, "setup_template": {"preset": "KS"}})

    base_id = CampaignGenerator().validate(base, base_dir=tmp_path).entry_ids["case"]
    ks_id = CampaignGenerator().validate(with_ks, base_dir=tmp_path).entry_ids["case"]

    assert with_ks.to_dict()["setup_template"]["preset"] == "KS"
    assert base_id != ks_id


def test_entry_metadata_remains_identity_bearing(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    first = Campaign.from_dict({
        "name": "metadata-identity",
        "entries": [{"name": "case", "source": "bcc", "metadata": {"recipe": "A"}}],
    })
    second = Campaign.from_dict({
        "name": "metadata-identity",
        "entries": [{"name": "case", "source": "bcc", "metadata": {"recipe": "B"}}],
    })

    first_id = CampaignGenerator().validate(first, base_dir=tmp_path).entry_ids["case"]
    second_id = CampaignGenerator().validate(second, base_dir=tmp_path).entry_ids["case"]

    assert first_id != second_id


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


@pytest.mark.parametrize("filename", ["lattice_new.nml", "input_new.nml"])
def test_reuse_rejects_modified_materialized_input(tmp_path: Path, filename: str) -> None:
    _source(tmp_path, "bccFe", "bcc")
    _source(tmp_path, "FeN", "fen")
    campaign = _campaign(tmp_path)
    generator = CampaignGenerator(setup_service=SetupService(runner=_FakeRunner()))
    destination = tmp_path / "output"
    generated = generator.generate(campaign, destination, base_dir=tmp_path)
    run_dir = destination / generated.entries[0].run_path
    path = run_dir / filename
    path.write_text(path.read_text(encoding="utf-8") + "\n! tampered\n", encoding="utf-8")

    with pytest.raises(CampaignCollisionError, match="INPUT_MODIFIED"):
        generator.generate(campaign, destination, base_dir=tmp_path)


def test_reuse_rejects_missing_materialized_input(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    _source(tmp_path, "FeN", "fen")
    campaign = _campaign(tmp_path)
    generator = CampaignGenerator(setup_service=SetupService(runner=_FakeRunner()))
    destination = tmp_path / "output"
    generated = generator.generate(campaign, destination, base_dir=tmp_path)
    run_dir = destination / generated.entries[0].run_path
    (run_dir / "input_new.nml").unlink()

    with pytest.raises(CampaignCollisionError, match="OWNED_FILE_MISSING"):
        generator.generate(campaign, destination, base_dir=tmp_path)


def test_runtime_outputs_do_not_change_materialized_input_identity(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    _source(tmp_path, "FeN", "fen")
    campaign = _campaign(tmp_path)
    generator = CampaignGenerator(setup_service=SetupService(runner=_FakeRunner()))
    destination = tmp_path / "output"
    generated = generator.generate(campaign, destination, base_dir=tmp_path)
    run_dir = destination / generated.entries[0].run_path
    (run_dir / "report.out").write_text("runtime output\n", encoding="utf-8")

    reused = generator.generate(campaign, destination, base_dir=tmp_path)

    assert reused.entries[0].entry_id == generated.entries[0].entry_id


def _valid_nested_exchange_manifest() -> dict:
    spec = {
        "pair_keys": [{"i": 0, "j": 0, "R": [0, 0, 0]}],
        "legacy_ijpair": [[1], [1]],
    }
    from rsutil.core.hashing import hash_json

    return {
        "schema_version": "rsutil-exchange-v1",
        "owner": "rsutil",
        "operation": "exchange_job",
        "spec": spec,
        "spec_hash": hash_json(spec),
        "canonical_pair_list": spec["pair_keys"],
        "pair_keys": spec["pair_keys"],
        "legacy_ijpair": spec["legacy_ijpair"],
    }


def test_parent_reuse_ignores_valid_nested_exchange_inputs(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    _source(tmp_path, "FeN", "fen")
    campaign = _campaign(tmp_path)
    generator = CampaignGenerator(setup_service=SetupService(runner=_FakeRunner()))
    destination = tmp_path / "output"
    generated = generator.generate(campaign, destination, base_dir=tmp_path)
    run_dir = destination / generated.entries[0].run_path
    child = run_dir / "Atom_1"
    child.mkdir()
    (child / "input.nml").write_text("&calculation\n post_processing='exchange'\n/\n", encoding="utf-8")
    (child / "rsutil_exchange_manifest.json").write_text(
        json.dumps(_valid_nested_exchange_manifest()), encoding="utf-8"
    )

    reused = generator.generate(campaign, destination, base_dir=tmp_path)
    assert reused.entries[0].entry_id == generated.entries[0].entry_id


def test_parent_reuse_rejects_unowned_nested_scientific_input(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    _source(tmp_path, "FeN", "fen")
    campaign = _campaign(tmp_path)
    generator = CampaignGenerator(setup_service=SetupService(runner=_FakeRunner()))
    destination = tmp_path / "output"
    generated = generator.generate(campaign, destination, base_dir=tmp_path)
    run_dir = destination / generated.entries[0].run_path
    child = run_dir / "unknown"
    child.mkdir()
    (child / "input.nml").write_text("&calculation\n/\n", encoding="utf-8")

    with pytest.raises(CampaignCollisionError, match="UNOWNED_DESCENDANT"):
        generator.generate(campaign, destination, base_dir=tmp_path)


def test_reuse_rejects_coordinated_input_and_historical_manifest_tampering(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    _source(tmp_path, "FeN", "fen")
    campaign = _campaign(tmp_path)
    generator = CampaignGenerator(setup_service=SetupService(runner=_FakeRunner()))
    destination = tmp_path / "output"
    generated = generator.generate(campaign, destination, base_dir=tmp_path)
    run_dir = destination / generated.entries[0].run_path
    manifest_path = destination / generated.entries[0].manifest_path
    tampered = run_dir / "input_new.nml"
    tampered.write_text("&atoms\n label(1) = 'Ni'\n/\n", encoding="utf-8")

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    record = {"sha256": hash_file(tampered), "size": tampered.stat().st_size}
    payload["output_hashes"]["input_new.nml"] = record
    payload["campaign"]["materialized_inputs"]["input_new.nml"] = record
    payload["campaign"]["integrity_hash"] = _campaign_integrity_hash(payload)
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignCollisionError, match="INPUT_MODIFIED"):
        generator.generate(campaign, destination, base_dir=tmp_path)


def test_seeded_stochastic_generation_can_be_reused(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    _source(tmp_path, "FeN", "fen")
    campaign = _campaign(tmp_path, second_setup={"lattice_dist": 0.01, "seed": 42})
    generator = CampaignGenerator(setup_service=SetupService(runner=_FakeRunner()))
    destination = tmp_path / "output"

    generated = generator.generate(campaign, destination, base_dir=tmp_path)
    reused = generator.generate(campaign, destination, base_dir=tmp_path)

    assert [entry.entry_id for entry in reused.entries] == [entry.entry_id for entry in generated.entries]


def test_unseeded_stochastic_campaign_is_rejected(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    _source(tmp_path, "FeN", "fen")
    campaign = _campaign(tmp_path, second_setup={"lattice_dist": 0.01})

    with pytest.raises(CampaignError, match="preflight validation failed"):
        CampaignGenerator(setup_service=SetupService(runner=_FakeRunner())).generate(
            campaign, tmp_path / "output", base_dir=tmp_path
        )


def test_reuse_rejects_tampered_campaign_manifest_metadata(tmp_path: Path) -> None:
    _source(tmp_path, "bccFe", "bcc")
    _source(tmp_path, "FeN", "fen")
    campaign = _campaign(tmp_path)
    generator = CampaignGenerator(setup_service=SetupService(runner=_FakeRunner()))
    destination = tmp_path / "output"
    generated = generator.generate(campaign, destination, base_dir=tmp_path)
    manifest_path = destination / generated.entries[0].manifest_path
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["campaign"]["calculation"]["mode"] = "tampered"
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignCollisionError, match="MANIFEST_MISMATCH"):
        generator.generate(campaign, destination, base_dir=tmp_path)
