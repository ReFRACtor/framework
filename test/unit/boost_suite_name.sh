#!/bin/bash

# Parses out test suite names from source files
# Outputs as a list seperated by ; so that CMake interprets as a list
grep '^[ ]*BOOST_FIXTURE_TEST_SUITE' $* | sed 's/.*(\([[:alnum:]_]*\),.*)/\1/' | tr '\n' ';'
