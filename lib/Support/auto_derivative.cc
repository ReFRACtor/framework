#include "auto_derivative.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class T> template<class Archive>
void AutoDerivative<T>::serialize(Archive & ar,
			      const unsigned int UNUSED(version))
{
  boost::serialization::void_cast_register<AutoDerivative<T>,
					   GenericObject>();
  ar & FP_NVP(val) & FP_NVP(grad);
}

FP_IMPLEMENT(AutoDerivativeDouble);
#endif

