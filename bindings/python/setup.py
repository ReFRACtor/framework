import os
from setuptools import setup, find_namespace_packages


def release_version():
    """The release number from version.py, which CMake writes next to this file from REFRACTOR_MAJOR_VERSION."""
    ns = {}
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "version.py")) as fh:
        exec(fh.read().split("SOURCE_DIRECTORY")[0], ns)
    return ns["__version__"]


# Namespace packages are a bit on the new side. If you haven't seen
# this before, look at https://packaging.python.org/guides/packaging-namespace-packages/#pkgutil-style-namespace-packages for a description of
# this.
setup(
    name='refractor-framework-swig',
    version=release_version(),
    description='ReFRACtor Framework SWIG code',
    author='James McDuffie',
    author_email='James.McDuffie@jpl.nasa.gov',
    packages=find_namespace_packages(include=["refractor.*"],),
    install_requires=[
        "refractor-framework",
    ],
    package_data={'' : ['_swig_wrap.so']},
)
