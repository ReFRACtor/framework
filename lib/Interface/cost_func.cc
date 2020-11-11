#include <cost_func.h>
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void CostFunc::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ProblemState)
    & FP_NVP(mssg) & FP_NVP(c_count);
}

FP_IMPLEMENT(CostFunc);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(CostFunc)
REGISTER_LUA_END()
#endif



const char * CostFunc::message_str() const
{
  switch( mssg ) {
  case SOLVED: return "SOLVED";
  case ERROR: return "ERROR";
  default: return "NONE";
  }
}
