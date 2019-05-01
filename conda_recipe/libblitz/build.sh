./configure --prefix="$PREFIX" --enable-shared --disable-doxygen
make
make install
rm -f "$PREFIX/lib/pkgconfig/blitz-uninstalled.pc"
rm -f "$PREFIX/lib/libblitz.la"

