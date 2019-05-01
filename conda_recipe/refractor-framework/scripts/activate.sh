#!/bin/bash

# Store existing env vars and set to this conda env

if [[ -n "$REFRACTOR_INPUTS" ]]; then
    export _CONDA_SET_REFRACTOR_INPUTS=$REFRACTOR_INPUTS
fi

if [[ -n "$LUA_PATH" ]]; then
    export _CONDA_SET_LUA_PATH=$LUA_PATH
fi

export REFRACTOR_INPUTS=${CONDA_PREFIX}/etc/refractor/input
if [ -n "$LUA_PATH" ]; then
    export LUA_PATH="${CONDA_PREFIX}/etc/refractor/config/?.lua;$LUA_PATH"
else
    export LUA_PATH="${CONDA_PREFIX}/etc/refractor/config/?.lua;"
fi

