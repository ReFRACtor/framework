#include "array_ad_with_unit.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class T, int D> template<class Archive>
void ArrayAdWithUnit<T, D>::serialize(Archive & ar,
					   const unsigned int UNUSED(version))
{
  boost::serialization::void_cast_register<ArrayAdWithUnit<T, D>,
					   GenericObject>();
  ar & FP_NVP(value) & FP_NVP(units);
}

FP_IMPLEMENT(ArrayAdWithUnit_double_1);
FP_IMPLEMENT(ArrayAdWithUnit_double_2);
FP_IMPLEMENT(ArrayAdWithUnit_double_3);
FP_IMPLEMENT(ArrayAdWithUnit_double_4);
#endif
