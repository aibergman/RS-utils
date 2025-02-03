import os
import subprocess

def build_wheel():
    """Custom build process for Poetry (runs CMake)"""
    print("🔧 Running CMake to build Fortran extension...")

    # Ensure the build directory exists
    os.makedirs("build", exist_ok=True)
    os.chdir("build")

    # Run CMake with Python bindings enabled
    subprocess.run(["cmake", "-DENABLE_PYTHON_BINDINGS=ON", ".."], check=True)
    subprocess.run(["make"], check=True)

    print("✅ CMake build complete. Compiled shared library should be in rsutil/_rsutil.so")

if __name__ == "__main__":
    build_wheel()

