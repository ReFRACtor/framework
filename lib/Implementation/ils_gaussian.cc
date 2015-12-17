#include "ils_gaussian.h"
using namespace FullPhysics;
using namespace blitz;

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
