import shutil
import pathlib
import subprocess

def copy_so_file():
    """Copy the compiled _rsutil.so file from build to rsutil/"""
    build_dir = pathlib.Path("build")
    target_dir = pathlib.Path(__file__).parent

    # Find the compiled shared library (matches _rsutil.*.so)
    so_files = list(build_dir.glob("_rsutil*.so"))

    if not so_files:
        raise FileNotFoundError("Could not find _rsutil.so in the build directory.")

    # Rename and copy the file
    shutil.copy(so_files[0], target_dir / "_rsutil.so")
    print(f"✅ Copied {so_files[0]} to {target_dir / '_rsutil.so'}")

def run_cmake():
    """Run CMake with Python bindings enabled"""
    print("🔧 Running CMake with Python bindings enabled...")
    subprocess.run(["cmake", "-DENABLE_PYTHON_BINDINGS=ON", "."], check=True)
    subprocess.run(["make"], check=True)
    copy_so_file()
