"""HT-R08 packaging and optional-import regression coverage."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_cli_core_import_does_not_require_visualization_dependency() -> None:
    script = """
import sys

class BlockPyVista:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'pyvista' or fullname.startswith('pyvista.'):
            raise ModuleNotFoundError('blocked for optional-dependency regression')
        return None

sys.meta_path.insert(0, BlockPyVista())
import rsutil.cli  # noqa: F401
assert 'pyvista' not in sys.modules
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"core CLI import unexpectedly needs PyVista\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
