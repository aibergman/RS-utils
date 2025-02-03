import os
import sys
import ctypes

# Dynamically determine the correct shared library filename
lib_filename = None
for filename in os.listdir(os.path.dirname(__file__)):
    print(filename)
    if filename.startswith("rsutil") and filename.endswith(".so"):
        lib_filename = filename
        break

if lib_filename is None:
    raise ImportError("Could not find the rsutil shared library.")

# Load the shared library
lib_path = os.path.join(os.path.dirname(__file__), lib_filename)
print('Loading shared library:', lib_path)
_rsutil = ctypes.CDLL(lib_path)

# Define process_lattice() signature
_rsutil.process_lattice.argtypes = []  # Modify if needed
_rsutil.process_lattice.restype = None  # Modify if needed

def process_lattice():
    """Python wrapper for the process_lattice() Fortran routine."""
    _rsutil.process_lattice()

