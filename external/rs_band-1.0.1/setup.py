#!/usr/bin/env python3
"""
Setup script for rs_band package
"""

from setuptools import setup, find_packages
import pathlib

# Read the README file
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

# Read requirements
requirements = []
try:
    with open("requirements.txt", "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                requirements.append(line)
except FileNotFoundError:
    # Fallback requirements
    requirements = [
        "numpy>=1.20.0",
        "matplotlib>=3.5.0", 
        "scipy>=1.7.0",
        "spglib>=2.0.0"
    ]

setup(
    name="rs_band",
    version="1.0.0",
    author="TB-LMTO Analysis Suite",
    author_email="contact@example.com",
    description="Real-space to Band structure analysis for TB-LMTO",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/username/rs_band",
    project_urls={
        "Bug Tracker": "https://github.com/username/rs_band/issues",
        "Documentation": "https://github.com/username/rs_band/blob/main/README.md",
        "Source Code": "https://github.com/username/rs_band",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    keywords="electronic-structure, band-structure, density-of-states, tb-lmto, physics, materials-science",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "full": ["pymatgen>=2022.0.0"],
        "dev": ["pytest>=6.0.0", "black", "flake8"],
        "docs": ["sphinx", "sphinx-rtd-theme"],
    },
    entry_points={
        'console_scripts': [
            'rs-band=rs_band.cli:main',
            'rsband=rs_band.cli:main',
            'rs-band-quick=rs_band.cli:quick_analysis',
            'rs-band-standard=rs_band.cli:standard_analysis', 
            'rs-band-high=rs_band.cli:high_analysis',
            'rs-band-structure=rs_band.cli:structure_analysis',
            'rs-band-dos-moments=rs_band.cli:dos_moments_analysis',
        ],
    },
    include_package_data=True,
    package_data={
        "rs_band": ["examples/*", "docs/*"],
    },
    zip_safe=False,
)
