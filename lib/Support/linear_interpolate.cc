#include "linear_interpolate.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class TX, class TY> template<class Archive>
void InterpolatePoint<TX, TY>::serialize(Archive & UNUSED(ar),
					 const unsigned int UNUSED(version))
{
    // Nothing to do
}

template<class TX, class TY> template<class Archive>
void Return1Point<TX, TY>::serialize(Archive & ar,
					     const unsigned int UNUSED(version))
{
  boost::serialization::void_cast_register<Return1Point<TX, TY>,
					   InterpolatePoint<TX, TY> >();
  ar & FP_NVP_(x0) & FP_NVP_(y0);
}

template<class TX, class TY> template<class Archive>
void LinearInterpolate2Point<TX, TY>::serialize(Archive & ar,
					     const unsigned int UNUSED(version))
{
  boost::serialization::void_cast_register<LinearInterpolate2Point<TX, TY>,
					   InterpolatePoint<TX, TY> >();
  ar & FP_NVP_(x0) & FP_NVP_(x1) & FP_NVP_(y0) & FP_NVP_(delta_y0);
}

template<class TX, class TY> template<class Archive>
void LinearInterpolate<TX, TY>::serialize(Archive & ar,
					     const unsigned int UNUSED(version))
{
  boost::serialization::void_cast_register<LinearInterpolate<TX, TY>,
					   GenericObject>();
  ar & FP_NVP(inter) & FP_NVP(out_of_range);
}

FP_IMPLEMENT(InterpolatePoint_double_double);
FP_IMPLEMENT(Return1Point_double_double);
FP_IMPLEMENT(LinearInterpolate2Point_double_double);
FP_IMPLEMENT(LinearInterpolate_double_double);


// See comments in linear_interpolate.h
// FP_IMPLEMENT(InterpolatePoint_double_auto_derivative_double);
// FP_IMPLEMENT(Return1Point_double_auto_derivative_double);
// FP_IMPLEMENT(LinearInterpolate2Point_double_auto_derivative_double);
// FP_IMPLEMENT(LinearInterpolate_double_auto_derivative_double);

#endif
