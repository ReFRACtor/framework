#!/usr/bin/env bash

if [ -d ~/refractor/absco ] 
then
    echo "Found ~/refractor/absco. Creating a link at ./absco" 
    ln -snf ~/refractor/absco absco
else
    echo "Directory ~/refractor/absco does not exists. This is not critical. However you won't be able to run / debug the Python tests."
fi
