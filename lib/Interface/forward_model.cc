#include "forward_model.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(ForwardModel)
.def("setup_grid", &ForwardModel::setup_grid)
REGISTER_LUA_END()
#endif
