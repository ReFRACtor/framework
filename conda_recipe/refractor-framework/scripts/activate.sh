#!/bin/bash

# Store existing env vars and set to this conda env

if [[ -n "$REFRACTOR_INPUT_PATH" ]]; then
    export _CONDA_SET_REFRACTOR_INPUT_PATH=$REFRACTOR_INPUT_PATH
fi

if [[ -n "$REFRACTOR_TEMPLATE_PATH" ]]; then
    export _CONDA_SET_REFRACTOR_TEMPLATE_PATH=$REFRACTOR_TEMPLATE_PATH
fi

if [[ -n "$LUA_PATH" ]]; then
    export _CONDA_SET_LUA_PATH=$LUA_PATH
fi

export REFRACTOR_INPUT_PATH=${CONDA_PREFIX}/etc/refractor/input
export REFRACTOR_TEMPLATE_PATH=${CONDA_PREFIX}/etc/refractor/template

if [ -n "$LUA_PATH" ]; then
    export LUA_PATH="${CONDA_PREFIX}/etc/refractor/config/?.lua;$LUA_PATH"
else
    export LUA_PATH="${CONDA_PREFIX}/etc/refractor/config/?.lua;"
fi

