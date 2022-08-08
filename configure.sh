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

pip install -r requirements.txt