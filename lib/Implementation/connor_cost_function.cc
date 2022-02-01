#include "connor_cost_function.h"

#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ConnorCostFunction, CostFunction)
.def(luabind::constructor<const boost::shared_ptr<StateVector>&, const boost::shared_ptr<ForwardModel>&, const boost::shared_ptr<Observation>&>())
REGISTER_LUA_END()
#endif

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ConnorCostFunction::serialize(Archive& ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CostFunction)
     & FP_NVP(statev) & FP_NVP(forward_model) & FP_NVP(meas);
}

FP_IMPLEMENT(ConnorCostFunction);
#endif
