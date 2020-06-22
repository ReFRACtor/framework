#include "array_ad.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<> template<class Archive>
void ArrayAd<double, 1>::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(ArrayAd_double_1);
  ar & FP_NVP(val) & FP_NVP(jac) & FP_NVP(is_const);
}

template<> template<class Archive>
void ArrayAd<double, 2>::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(ArrayAd_double_2);
  ar & FP_NVP(val) & FP_NVP(jac) & FP_NVP(is_const);
}

template<> template<class Archive>
void ArrayAd<double, 3>::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(ArrayAd_double_3);
  ar & FP_NVP(val) & FP_NVP(jac) & FP_NVP(is_const);
}

template<> template<class Archive>
void ArrayAd<double, 4>::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(ArrayAd_double_4);
  ar & FP_NVP(val) & FP_NVP(jac) & FP_NVP(is_const);
}

FP_IMPLEMENT(ArrayAd_double_1);
FP_IMPLEMENT(ArrayAd_double_2);
FP_IMPLEMENT(ArrayAd_double_3);
FP_IMPLEMENT(ArrayAd_double_4);
#endif
