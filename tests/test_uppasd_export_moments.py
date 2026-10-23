"""UppASD momfile moments come from report.out, with a warned fallback."""

from __future__ import annotations

from pathlib import Path

import pytest

from rsutil.export.uppasd import export_from_files, read_report_spin_moments

LATTICE = """&lattice
nbulk_bulk = 2
ntot = 2
nbas = 2
nrec = 2
a(:, 1) = 1.0, 0.0, 0.0
a(:, 2) = 0.0, 1.0, 0.0
a(:, 3) = 0.0, 0.0, 1.0
crd(:, 1) = 0.0, 0.0, 0.0
crd(:, 2) = 0.5, 0.5, 0.5
/
"""

REPORT = """Total spin moment:    -0.0 -0.0 4.5
Spin moment of atom   1:  2.778656
Spin moment projections of atom   1: -0.000000 -0.000000  2.778656
Spin moment of atom   2:  1.760090
Spin moment projections of atom   2: -0.000000 -0.000000  1.760090
Orbital moment of atom   1:  0.052501
"""


def _moments(out_dir: Path) -> list[float]:
    return [float(line.split()[2]) for line in (out_dir / "momfile").read_text().splitlines()]


def _setup(tmp_path: Path, report: str | None) -> tuple[Path, Path | None]:
    lattice = tmp_path / "lattice.nml"
    lattice.write_text(LATTICE)
    report_path = None
    if report is not None:
        report_path = tmp_path / "report.out"
        report_path.write_text(report)
    return lattice, report_path


def test_read_report_spin_moments_ignores_orbital_and_projection_lines(tmp_path: Path) -> None:
    _, report = _setup(tmp_path, REPORT)
    assert read_report_spin_moments(str(report)) == {1: 2.778656, 2: 1.76009}


def test_momfile_uses_report_moments(tmp_path: Path) -> None:
    lattice, report = _setup(tmp_path, REPORT)
    out = tmp_path / "out"
    export_from_files(lattice_nml=str(lattice), out_dir=str(out), report_file=str(report))
    assert _moments(out) == [2.778656, 1.76009]


def test_missing_report_falls_back_to_default_with_warning(tmp_path: Path) -> None:
    lattice, _ = _setup(tmp_path, None)
    out = tmp_path / "out"
    with pytest.warns(UserWarning, match="report.out not found"):
        export_from_files(
            lattice_nml=str(lattice), out_dir=str(out), report_file=str(tmp_path / "report.out")
        )
    assert _moments(out) == [1.0, 1.0]
    with pytest.warns(UserWarning, match="report.out not found"):
        export_from_files(lattice_nml=str(lattice), out_dir=str(out))
    assert _moments(out) == [1.0, 1.0]


def test_atom_missing_from_report_falls_back_individually(tmp_path: Path) -> None:
    lattice, report = _setup(tmp_path, "Spin moment of atom   1:  2.5\n")
    out = tmp_path / "out"
    with pytest.warns(UserWarning, match=r"atom\(s\) \[2\]"):
        export_from_files(lattice_nml=str(lattice), out_dir=str(out), report_file=str(report))
    assert _moments(out) == [2.5, 1.0]


def test_cli_defaults_report_to_data_directory(tmp_path: Path, capsys) -> None:
    import sys
    from rsutil.cli import main

    _setup(tmp_path, REPORT)
    old = sys.argv
    sys.argv = ["rsutil", "--export", "--data", str(tmp_path)]
    try:
        main()
    finally:
        sys.argv = old
    assert _moments(tmp_path / "uppasd_out") == [2.778656, 1.76009]
