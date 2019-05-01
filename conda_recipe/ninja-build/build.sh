
#!/bin/bash

mkdir -p "${PREFIX}"
mkdir -p "${PREFIX}/bin"
mkdir -p "${PREFIX}/include"
mkdir -p "${PREFIX}/lib"

./configure.py --bootstrap

cp -p ninja "$PREFIX/bin/ninja-build"
