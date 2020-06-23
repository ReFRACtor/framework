#include "array_ad.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class T, int D> template<class Archive>
void ArrayAd<T, D>::serialize(Archive & ar,
			      const unsigned int UNUSED(version))
{
  boost::serialization::void_cast_register<ArrayAd<T, D>,
					   GenericObject>();
  ar & FP_NVP(val) & FP_NVP(jac) & FP_NVP(is_const);
}

FP_IMPLEMENT(ArrayAd_double_1);
FP_IMPLEMENT(ArrayAd_double_2);
FP_IMPLEMENT(ArrayAd_double_3);
FP_IMPLEMENT(ArrayAd_double_4);
#endif
