#include "scattering_moment_interpolator.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ScatteringMomentInterpolate2Point::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  boost::serialization::void_cast_register<ScatteringMomentInterpolate2Point,
					   GenericObject>();
  ar & FP_NVP(wn0) & FP_NVP(pf0) & FP_NVP(delta_pf0);
}

template<class Archive>
void ScatteringMomentInterpolate::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  boost::serialization::void_cast_register<ScatteringMomentInterpolate,
					   GenericObject>();
  ar & FP_NVP(inter);
}

FP_IMPLEMENT(ScatteringMomentInterpolate2Point);
FP_IMPLEMENT(ScatteringMomentInterpolate);
#endif
