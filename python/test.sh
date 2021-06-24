#! /bin/bash

PYTHONPATH=$(pwd):${PYTHONPATH} py.test -rxXs --log-cli-level=info test
