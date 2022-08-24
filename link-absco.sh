#!/usr/bin/env bash

if [ -d ~/refractor-dev/absco ] 
then
    echo "Found ~/refractor-dev/absco. Creating a link at ./absco" 
    ln -snf ~/refractor-dev/absco absco
else
    if [ -d ~/refractor/absco ] 
    then
        echo "Found ~/refractor/absco ... Creating a link at ./absco" 
        ln -snf ~/refractor/absco absco
    else
        echo "Directories ~/refractor-dev/absco or ~/refractor/absco do not exist. \
              This is not critical. However you won't be able to run / debug the Python tests."
    fi

fi
