#include "forward_model.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ForwardModel::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StackedRadianceMixin);
}

FP_IMPLEMENT(ForwardModel);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(ForwardModel)
.def("setup_grid", &ForwardModel::setup_grid)
REGISTER_LUA_END()
#endif
