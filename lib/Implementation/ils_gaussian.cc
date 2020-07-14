#include "ils_gaussian.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void IlsGaussian::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IlsFunction)
    & FP_NVP(gamma) & FP_NVP_(band_name) & FP_NVP_(hdf_band_name);
}

FP_IMPLEMENT(IlsGaussian);
#endif

void IlsGaussian::ils
(const AutoDerivative<double>& wn_center,
 const blitz::Array<double, 1>& wn,
 ArrayAd<double, 1>& res_a) const
{
  const double sqrt_PI = 1.7724538509055159;
  Array<AutoDerivative<double>, 1> res(wn.shape());
  res = exp(-sqr((wn - wn_center) / gamma)) / (gamma*sqrt_PI);
  res_a.reference(ArrayAd<double, 1>(res));
}
