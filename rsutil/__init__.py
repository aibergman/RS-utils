import sys
import importlib.util
import pathlib
import glob

# Dynamically locate the shared library file
module_name = "_rsutil"
module_path = pathlib.Path(__file__).parent

# Find the actual shared library file (supports different OS and architectures)
shared_libs = list(module_path.glob(f"{module_name}.cpython-{sys.version_info.major}{sys.version_info.minor}*.so"))

if not shared_libs:
    raise ImportError(f"Cannot find compiled shared library for {module_name} in {module_path}")

# Load the shared object file
spec = importlib.util.spec_from_file_location(module_name, shared_libs[0])
_rsutil = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_rsutil)

# 🔹 Fix the function name (since `f90wrap` renamed it)
process_lattice = _rsutil.f90wrap_process_lattice

