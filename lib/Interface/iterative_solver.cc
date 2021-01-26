#include <iterative_solver.h>
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void IterativeSolver::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableIterativeSolver)
    & FP_NVP(max_cost_f_calls) & FP_NVP(verbose) & FP_NVP(stat)
    & FP_NVP(Accepted_points) & FP_NVP(Cost_at_accepted_points);
}

FP_IMPLEMENT(IterativeSolver);
FP_OBSERVER_SERIALIZE(IterativeSolver);
#endif


#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(IterativeSolver)
REGISTER_LUA_END()
#endif



const char * IterativeSolver::status_str() const
{
  switch( stat ) {
  case SUCCESS: return "SUCCESS";
  case CONTINUE: return "CONTINUE";
  case STALLED: return "STALLED";
  case ERROR: return "ERROR";
  case UNTRIED: return "UNTRIED";
  }
  return "UNKNOWN";
}
