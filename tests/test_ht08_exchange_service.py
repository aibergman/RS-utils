"""HT-08 exchange-service contracts and failure propagation."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from rsutil.exchange.exchange_map import (
    find_matching_indices,
    get_neighbors_list,
    read_rs_clust,
)
from rsutil.exchange.fermi import FermiEnergyParseError, parse_fermi_energy
from rsutil.exchange.service import ExchangeService

from .conftest import copy_example_files, run_rsutil


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


def test_exchange_service_is_path_explicit_and_pair_inventory_matches_legacy_set(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source"
    source.mkdir()
    _copy_exchange_fixture(source)
    unrelated_cwd = tmp_path / "unrelated"
    unrelated_cwd.mkdir()
    monkeypatch.chdir(unrelated_cwd)

    service = ExchangeService(source)
    plan = service.plan(cutoff=1.0)

    cluster = read_rs_clust(source / "clust", source / "input.nml")
    legacy_neighbors = get_neighbors_list(cluster, plan.reference.recursion_sites, cutoff=1.0)
    assert [job.matching_indices for job in plan.jobs] == [
        tuple(indices)
        for indices in find_matching_indices(cluster, legacy_neighbors).values()
    ]
    assert plan.pairs
    assert all(type(pair.key).__name__ == "PairKey" for pair in plan.pairs)
    assert {pair.key for pair in plan.pairs} == set(plan.pair_keys)

    result = service.generate(cutoff=1.0, rij_out="mapping.json")
    assert result.mapping_path == source / "mapping.json"
    assert result.mapping_path.is_file()
    mapping = json.loads(result.mapping_path.read_text(encoding="utf-8"))
    assert mapping[0]["i"] in (1, 2)
    assert mapping[0]["j"] in plan.jobs[0].matching_indices + plan.jobs[1].matching_indices
    assert (source / "Atom_1" / "input.nml").is_file()
    assert (source / "Atom_2" / "input.nml").is_file()


@pytest.mark.parametrize(
    "contents",
    ["no Fermi record\n", "Fermi energy: not-a-number\n"],
)
def test_fermi_parser_returns_structured_error_and_cli_fails(
    tmp_path: Path, contents: str
) -> None:
    source = tmp_path / "source"
    source.mkdir()
    _copy_exchange_fixture(source)
    (source / "report.out").write_text(contents, encoding="utf-8")

    with pytest.raises(FermiEnergyParseError):
        parse_fermi_energy(source / "report.out")
    result = run_rsutil(source, "-J", "--cutoff", "1.0")
    assert result.returncode != 0
    assert "Exchange generation failed" in result.stderr
    assert not (source / "Atom_1").exists()


def test_explicit_rij_output_failure_fails_operation(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)
    blocked = tmp_path / "blocked"
    blocked.write_text("not a directory", encoding="utf-8")

    result = run_rsutil(
        tmp_path,
        "-J",
        "--cutoff",
        "1.0",
        "--rij-out",
        str(blocked / "mapping.json"),
    )
    assert result.returncode != 0
    assert "Exchange generation failed" in result.stderr
    assert not (tmp_path / "Atom_1").exists()


def test_full_pairs_opt_out_reproduces_unreduced_inventory(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)

    reduced = ExchangeService(tmp_path).plan(
        cutoff=1.0,
        unique=True,
        magnetic_reference="collinear_fm",
        magnetic_axis=(0.0, 0.0, 1.0),
    )
    full = ExchangeService(tmp_path).plan(
        cutoff=1.0, unique=True, full_pairs=True
    )

    assert len(full.pairs) > len(reduced.pairs)
    assert all(len(job.pairs) > 1 for job in full.jobs)


def test_unproven_afm_reference_disables_unique_reduction(tmp_path: Path) -> None:
    _copy_exchange_fixture(tmp_path)
    service = ExchangeService(tmp_path)

    afm = service.plan(cutoff=1.0, unique=True, magnetic_reference="afm")
    full = service.plan(cutoff=1.0, unique=True, full_pairs=True)

    assert set(afm.pair_keys) == set(full.pair_keys)
    assert all(job.magnetic_reference == "afm" for job in afm.jobs)
