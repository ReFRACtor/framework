import os
from setuptools import setup, find_namespace_packages

# Namespace packages are a bit on the new side. If you haven't seen
# this before, look at https://packaging.python.org/guides/packaging-namespace-packages/#pkgutil-style-namespace-packages for a description of
# this.
setup(
    name='refractor-framework-swig',
    version='1.0.0',
    description='ReFRACtor Framework SWIG code',
    author='James McDuffie',
    author_email='James.McDuffie@jpl.nasa.gov',
    packages=find_namespace_packages(include=["refractor.*"],),
    install_requires=[
        'numpy', 
    ],
    package_data={'' : ['_swig_wrap.so']},
)
