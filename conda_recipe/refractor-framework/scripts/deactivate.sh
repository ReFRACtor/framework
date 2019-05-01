#!/bin/bash
# Restore previous env vars if they were set

unset REFRACTOR_INPUTS
if [[ -n "$_CONDA_SET_REFRACTOR_INPUTS" ]]; then
    export REFRACTOR_INPUTS=$_CONDA_SET_REFRACTOR_INPUTS
    unset _CONDA_SET_REFRACTOR_INPUTS
fi

unset LUA_PATH
if [[ -n "$_CONDA_SET_LUA_PATH" ]]; then
    export LUA_PATH=$_CONDA_SET_LUA_PATH
    unset _CONDA_SET_LUA_PATH
fi

