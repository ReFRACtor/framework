import os
import re
from setuptools import setup, find_namespace_packages


def cmake_version():
    """The release number, REFRACTOR_MAJOR_VERSION in CMakeLists.txt, which is the one place it is set."""
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "CMakeLists.txt")) as fh:
        return re.search(r'set\(REFRACTOR_MAJOR_VERSION\s+"([^"]+)"\)', fh.read()).group(1)


# Namespace packages are a bit on the new side. If you haven't seen
# this before, look at https://packaging.python.org/guides/packaging-namespace-packages/#pkgutil-style-namespace-packages for a description of
# this.
setup(
    name='refractor-framework',
    version=cmake_version(),
    description='ReFRACtor Framework',
    author='James McDuffie',
    author_email='James.McDuffie@jpl.nasa.gov',
    packages=find_namespace_packages(include=["refractor.*"],
                                     where="python"),
    package_dir={"": "python"},
    install_requires=[
        "filelock",
        "future",
        "h5py",
        "jsonpickle",
        "matplotlib",
        "netcdf4",
        "nlopt",
        "numpy",
        "progressbar33",
        "pytest",
        "scipy",
        "setuptools",
        "six",
    ],
)
