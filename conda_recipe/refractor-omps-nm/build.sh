# Note that the various libraries are installed in $PREFIX which isn't
# searched, so point to these
cmake -DBUILD_TESTS=OFF -DCMAKE_INSTALL_PREFIX=$PREFIX $SRC_DIR -DCMAKE_INSTALL_LIBDIR=lib -DBOOST_ROOT=$PREFIX -DGSL_ROOT_DIR=$PREFIX -DBLITZ_DIR=$PREFIX -DFIRST_ORDER_DIR=$PREFIX -DLIDORT_DIR=$PREFIX -DLUABIND_ROOT_DIR=$PREFIX -DTWOSTREAM_DIR=$PREFIX -DREFRACTOR_DIR=$PREFIX
make -j${CPU_COUNT} all
make install

cat "$RECIPE_DIR"/extra_for_setup "$PREFIX"/setup_omps_nm_env.sh > "$PREFIX"/setup_omps_nm_env.sh.temp
mv -f "$PREFIX"/setup_omps_nm_env.sh.temp "$PREFIX"/setup_omps_nm_env.sh



