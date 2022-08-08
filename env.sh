#!/usr/bin/env bash

source configure.sh

# install our package for easy development
pip install --editable build/debug/bindings/python
pip install --editable .

# Inputs

# Link to home dir ~/refractor/absco
./link-absco.sh

export REFRACTOR_INPUT_PATH=$(realpath ./input)
export abscodir=$(realpath ./absco/unit_test_absco/tables)

# Generate .env file for pytest
echo "REFRACTOR_INPUT_PATH=$(realpath './input')" > .env
echo "abscodir=$(realpath './absco/unit_test_absco/tables')" >> .env
