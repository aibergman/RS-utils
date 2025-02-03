import sys
import importlib.util
import pathlib

# Dynamically load `_rsutil.so`
module_name = "_rsutil"
module_path = pathlib.Path(__file__).parent / f"{module_name}{sys.implementation.cache_tag}.so"

spec = importlib.util.spec_from_file_location(module_name, module_path)
_rsutil = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_rsutil)

process_lattice = _rsutil.process_lattice

