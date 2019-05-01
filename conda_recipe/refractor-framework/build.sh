# Note that the various libraries are installed in $PREFIX which isn't
# searched, so point to these
cmake -DBUILD_TESTS=OFF -DCMAKE_INSTALL_PREFIX=$PREFIX $SRC_DIR -DCMAKE_INSTALL_LIBDIR=lib -DBOOST_ROOT=$PREFIX -DGSL_ROOT_DIR=$PREFIX -DBLITZ_DIR=$PREFIX -DFIRST_ORDER_DIR=$PREFIX -DLIDORT_DIR=$PREFIX -DLUABIND_ROOT_DIR=$PREFIX -DTWOSTREAM_DIR=$PREFIX
make -j${CPU_COUNT} all
make install

# Replace any library in the build directory with the prefix directory
if [ `uname` = 'Linux' ]; then
    sed "s|${BUILD_PREFIX}|${PREFIX}|g" -i "$PREFIX/cmake/RefractorTargets.cmake";
else    
    sed -i '' "s|${BUILD_PREFIX}|${PREFIX}|g" "$PREFIX/cmake/RefractorTargets.cmake";
fi

cp "$RECIPE_DIR"/framework-post-link.sh "$PREFIX"/bin/.framework-post-link.sh
cp "$RECIPE_DIR"/framework-post-link.py "$PREFIX"/bin/.framework-post-link.py

cat "$RECIPE_DIR"/extra_for_setup "$PREFIX"/setup_fp_env.sh > "$PREFIX"/setup_fp_env.sh.temp
mv -f "$PREFIX"/setup_fp_env.sh.temp "$PREFIX"/setup_fp_env.sh

ACTIVATE_DIR=$PREFIX/etc/conda/activate.d
DEACTIVATE_DIR=$PREFIX/etc/conda/deactivate.d
mkdir -p $ACTIVATE_DIR
mkdir -p $DEACTIVATE_DIR

cp $RECIPE_DIR/scripts/activate.sh $ACTIVATE_DIR/framework-activate.sh
cp $RECIPE_DIR/scripts/deactivate.sh $DEACTIVATE_DIR/framework-deactivate.sh


