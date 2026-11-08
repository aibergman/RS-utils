# Reproducible HT testing

The supported test environment is Python 3.11 or newer. CI exercises the
complete Python 3.11 suite on Linux, from both an editable install and an
installed wheel.
The HT test extra contains pytest, ASE, and the optional
band-analysis dependencies used by the complete checked-in suite; the core
runtime does not require ASE or PyVista.

## Linux clean checkout

Install the system compiler and CMake first. Debian/Ubuntu example:

```bash
sudo apt-get update
sudo apt-get install -y gfortran cmake
```

Create an isolated environment and install the package with its test extra:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[test]'
python -m pytest -q
```

The editable install invokes the PEP 517 build backend and builds the
protected Fortran extension in the environment, placing the ABI-qualified
artifact directly in the rsutil/ package. The install itself is sufficient;
do not manually copy a shared library or run build_ext --inplace. To build a
wheel as a separate check, install the wheel builder and run:

```bash
python -m pip install build
python -m build --wheel
```

## Installed-wheel artifact verification

The artifact lane must install the built wheel into a venv outside the
checkout, with no editable install and no `PYTHONPATH` pointing at the source.
The following is the complete local verification sequence:

```bash
sudo apt-get update
sudo apt-get install -y gfortran cmake
python3 -m pip install --upgrade build
rm -rf wheelhouse
mkdir wheelhouse
python3 -m build --wheel --outdir wheelhouse
test "$(find wheelhouse -maxdepth 1 -type f -name '*.whl' | wc -l | tr -d ' ')" -eq 1

wheel_path="$(find wheelhouse -maxdepth 1 -type f -name '*.whl' -print -quit)"
wheel_venv="$(mktemp -d /tmp/rsutil-wheel-venv.XXXXXX)"
python3 -m venv "$wheel_venv"
"$wheel_venv/bin/python" -m pip install --upgrade pip
"$wheel_venv/bin/python" -m pip install "$PWD/$wheel_path[test]"
"$wheel_venv/bin/python" -m pip check

wheel_suite="$(mktemp -d /tmp/rsutil-wheel-suite.XXXXXX)"
cp -R tests examples "$wheel_suite/"
cd "$wheel_suite"
env -u PYTHONPATH PATH="$wheel_venv/bin:$PATH" "$wheel_venv/bin/python" - <<'PY'
import importlib.util
import pathlib
import site
import sys
import rsutil

site_roots = [pathlib.Path(path).resolve() for path in site.getsitepackages()]
package_path = pathlib.Path(rsutil.__file__).resolve()
native_spec = importlib.util.find_spec("rsutil.rsutil")
assert native_spec is not None and native_spec.origin is not None
native_path = pathlib.Path(native_spec.origin).resolve()
assert any(root == package_path or root in package_path.parents for root in site_roots)
assert any(root == native_path or root in native_path.parents for root in site_roots)
assert native_path.suffix in {".so", ".dylib", ".pyd"}
print(f"Python: {sys.executable}")
print(f"rsutil.__file__: {package_path}")
print(f"native extension: {native_path}")
print("site-packages assertion: PASS")
PY
env -u PYTHONPATH PATH="$wheel_venv/bin:$PATH" "$wheel_venv/bin/python" -m pytest -q "$wheel_suite/tests"
```

The external copy is intentional: it prevents checkout files from satisfying
imports or fixture lookups while the path assertion proves both the Python
package and compiled extension came from the installed wheel.

The release wheel therefore has the same package-local layout as the editable
install: the native artifact is embedded under rsutil/ and is loaded from
that installed package directory. The supported primary native release lane
is Linux; other platforms remain subject to their existing compiler and CMake
toolchain support.

The optional extras are:

- `ase` for `Atoms <-> RSStructure` interoperability;
- `bands` for band/DOS analysis and plotting;
- `visualization` for the PyVista-based structure viewer;
- `test` for the complete checked-in pytest suite;
- `full` for all optional runtime features.

The release candidate package version is `0.11.0rc1`. The wheel and installed
artifact checks in this page are the release evidence path; the editable
installation above remains the developer/source-test workflow.

## Setup fixed-point regression

The deterministic setup contract is tested by promoting `lattice_new.nml` and
`input_new.nml` from a first `-S` pass to the next pass with identical options.
The scientifically relevant generated namelists must be semantically equal;
formatting differences alone are not significant.  This regression records
the observed fixed-point property for representative single-site and
multisite/multispecies fixtures.  It is test evidence, not a user workflow:
production setup performs one pass only.
