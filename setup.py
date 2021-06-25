import os
from setuptools import setup, find_namespace_packages

# Namespace packages are a bit on the new side. If you haven't seen
# this before, look at https://packaging.python.org/guides/packaging-namespace-packages/#pkgutil-style-namespace-packages for a description of
# this.
setup(
    name='refractor-framework',
    version='1.0.0',
    description='ReFRACtor Framework',
    author='James McDuffie',
    author_email='James.McDuffie@jpl.nasa.gov',
    packages=find_namespace_packages(include=["refractor.*"],
                                     where="python"),
    package_dir={"": "python"},
    install_requires=[
        'numpy', 
    ],
)
