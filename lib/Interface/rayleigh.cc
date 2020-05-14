#include "rayleigh.h"

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Rayleigh)
.def("optical_depth_each_layer", &Rayleigh::optical_depth_each_layer)
REGISTER_LUA_END()
#endif
