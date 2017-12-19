#include "connor_cost_function.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ConnorCostFunction, CostFunction)
.def(luabind::constructor<const boost::shared_ptr<StateVector>&, const boost::shared_ptr<ForwardModel>&, const boost::shared_ptr<Observation>&>())
REGISTER_LUA_END()
#endif
