#include "spectral_range.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SpectralRange::serialize(Archive& ar,
			      const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(SpectralRange);
  ar & FP_NVP_(data) & FP_NVP_(units) & FP_NVP_(uncertainty);
}

FP_IMPLEMENT(SpectralRange);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
typedef const blitz::Array<double, 1>& (SpectralRange::*ftype)() const;
REGISTER_LUA_CLASS(SpectralRange)
.def("data", (ftype) &SpectralRange::data)
.def("units", &SpectralRange::units)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Convert to given units.
//-----------------------------------------------------------------------

SpectralRange SpectralRange::convert(const Unit& R) const
{
  double conv = FullPhysics::conversion(units_, R);
  Array<double, 1> dv(data_.value().shape()), uncer(uncertainty_.shape());
  Array<double, 2> djac(data_.jacobian().shape());
  dv = data_.value() * conv;
  if(djac.cols() > 0)
    djac = data_.jacobian() * conv;
  if(uncer.rows() > 0)
    uncer = uncertainty_ * conv;
  return SpectralRange(ArrayAd<double, 1>(dv, djac), R, uncer);
}
