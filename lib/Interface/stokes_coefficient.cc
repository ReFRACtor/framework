#include "stokes_coefficient.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StokesCoefficient::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableStokesCoefficient);
    
}

FP_IMPLEMENT(StokesCoefficient);
FP_OBSERVER_SERIALIZE(StokesCoefficient);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
blitz::Array<double, 2> stokes_coefficient_value(const StokesCoefficient& S)
{
  return S.stokes_coefficient().value();
}
REGISTER_LUA_CLASS(StokesCoefficient)
.def("stokes_coefficient", &StokesCoefficient::stokes_coefficient)
.def("stokes_coefficient_value", &stokes_coefficient_value)
REGISTER_LUA_END()
#endif

