#! /bin/bash

src_test_dir=$(dirname $0)/test

PYTHONPATH=$(pwd):${PYTHONPATH} py.test -rxXs --log-cli-level=info $src_test_dir
