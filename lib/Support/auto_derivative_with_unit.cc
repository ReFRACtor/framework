#include "auto_derivative_with_unit.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class T> template<class Archive>
void AutoDerivativeWithUnit<T>::serialize(Archive & ar,
			      const unsigned int UNUSED(version))
{
  boost::serialization::void_cast_register<AutoDerivativeWithUnit<T>,
					   GenericObject>();
  ar & FP_NVP(value) & FP_NVP(units);
}

FP_IMPLEMENT(AutoDerivativeWithUnitDouble);
#endif

