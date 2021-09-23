#include "pressure.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Pressure::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservablePressure);
    
}

FP_IMPLEMENT(Pressure);
FP_OBSERVER_SERIALIZE(Pressure);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Pressure)
.def("max_number_level", &Pressure::max_number_level)
.def("pressure_grid", &Pressure::pressure_grid)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Return surface pressure, which is just the pressure at the bottom
/// level of pressure_grid.
///
/// This is in Pascals.
//-----------------------------------------------------------------------

AutoDerivativeWithUnit<double> Pressure::surface_pressure() const
{
  ArrayAdWithUnit<double, 1> p(pressure_grid(INCREASING_PRESSURE));
  return p(p.rows() - 1);
}
