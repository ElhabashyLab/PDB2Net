# setup.py
from setuptools import setup, find_packages

setup(
    name="your_package",
    version="0.1",
    description="A package for reading .pdb, .cif, .mmcif, and .csv files.",
    packages=find_packages(),
    python_requires=">=3.6",
)
