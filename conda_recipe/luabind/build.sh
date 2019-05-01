# For GCC 7.3, luabind has a segmentation fault when compiled optimized. This
# is probably some bug in luabind that we could fix, but 1) luabind isn't
# really updated anymore and 2) we will move lua out of ReFactor in favor of
# python. So it isn't worth trying to fix. Just turn off optimization for newer
# versions of the compiler. Since this is just used by lua configuration, this
# shouldn't have much of an impact on over all performance
export CXXFLAGS="${CXXFLAGS} -O0"
cmake -DBUILD_TESTS=OFF -DCMAKE_INSTALL_PREFIX=$PREFIX $SRC_DIR -DCMAKE_INSTALL_LIBDIR=lib -DBOOST_ROOT="$CONDA_PREFIX"
make VERBOSE=1 luabind


