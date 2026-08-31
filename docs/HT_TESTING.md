# Reproducible HT testing

The supported test environment is Python 3.9 or newer. CI currently exercises
Python 3.11 on Linux. The HT test extra contains pytest, ASE, and the optional
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
protected Fortran extension in the environment. To build a wheel as a
separate check, install the wheel builder and run:

```bash
python -m pip install build
python -m build --wheel
```

The optional extras are:

- `ase` for `Atoms <-> RSStructure` interoperability;
- `bands` for band/DOS analysis and plotting;
- `visualization` for the PyVista-based structure viewer;
- `test` for the complete checked-in pytest suite;
- `full` for all optional runtime features.

The package remains at version `0.10.0` for this remediation. No release is
being cut by R08, so a version bump is deferred to the next release decision.
