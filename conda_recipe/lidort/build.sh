# Note there isn't any agreement I could find on where fortran modules
# should go. There is a lot of disagreement on the web, and some suggestions
# that they shouldn't even be used but rather have interface blocks. But
# we are stuck with the code we have. Install these in include for now,
# we can revisit this if it becomes an issue

# Needed if we clone from git, but currently we can't do that because
# conda doesn't support lfs
# SRC_DIR is actually a subdirectory of the full framework repository
#SRC_DIR="$SRC_DIR/thirdparty/lidort-3.7"
cmake -DBUILD_TESTS=OFF -DCMAKE_INSTALL_PREFIX=$PREFIX $SRC_DIR -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_Fortran_MODULE_DIRECTORY=$PREFIX/include
make -j${CPU_COUNT} all
make install
cd $PREFIX/lib

