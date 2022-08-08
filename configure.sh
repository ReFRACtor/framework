#!/usr/bin/env bash

# Python

# Do not install Python if in conda environment or we are on Docker container
if [[ -z "$CONDA_DEFAULT_ENV" && ! -f /.dockerenv ]]
then
    # install Python 3.8.0 if not installed
    pyenv install 3.8.0 --skip-existing
    pyenv versions

    # use Python 3 from .python-version for local development
    eval "$(pyenv init -)"
fi

# create virtual environment
if [ ! -z "$VIRTUAL_ENV" ]; then
    deactivate
fi
python3 -m venv .venv

# activate virtual environment
source .venv/bin/activate

# upgrade pip
pip install --upgrade pip wheel

# install runtime packages
pip install -r requirements.txt

# TODO: Make this automatic. Right now requires `make install` to be run before. 
# install our package for easy development
# pip install --editable build/bindings/python
# pip install --editable .

# Inputs

# Link to home dir ~/refractor/absco
./link-absco.sh

echo "REFRACTOR_INPUT_PATH=$(realpath './input')" > .env
echo "abscodir=$(realpath './absco/unit_test_absco/tables')" >> .env
