#include "stokes_coefficient_fraction.h"
#include "fp_serialize_support.h"
#include "ostream_pad.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StokesCoefficientFraction::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StokesCoefficientImpBase)
    & FP_NVP(stokes_coeff_parallel);
}

FP_IMPLEMENT(StokesCoefficientFraction);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(StokesCoefficientFraction, StokesCoefficient)
.def(luabind::constructor<const blitz::Array<double, 2>&, 
                          const blitz::Array<double, 1>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
StokesCoefficientFraction::StokesCoefficientFraction
(const blitz::Array<double, 2>& Stokes_coeff_parallel,
 const blitz::Array<double, 1>& Coeffs)
: stokes_coeff_parallel(Stokes_coeff_parallel.copy())
{
  if(stokes_coeff_parallel.rows() != Coeffs.rows())
    throw Exception("Stokes_coeff_parallel and Coeff need to have the same number of rows");
  stokes_coeff.resize(stokes_coeff_parallel.shape(), 0);
  stokes_coeff.value() = stokes_coeff_parallel;
  for(int i = 0; i < stokes_coeff.rows(); ++i) {
    stokes_coeff.value()(i, 1) *= (1 - 2 * Coeffs(i));
    stokes_coeff.value()(i, 2) *= (1 - 2 * Coeffs(i));
  }
  init(Coeffs);
}

void StokesCoefficientFraction::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "StokesCoefficientFraction:\n";
  Os << "  Initial parallel stokes coefficients:\n";
  opad << stokes_coeff_parallel << "\n";
  opad.strict_sync();
  Os << "  Coefficient:\n";
  opad << coeff.value() << "\n";
  opad.strict_sync();
  Os << "  Stokes coefficient:\n";
  opad << stokes_coefficient().value();
  opad.strict_sync();
}

boost::shared_ptr<StokesCoefficient> StokesCoefficientFraction::clone() const
{
  return boost::shared_ptr<StokesCoefficient>
    (new StokesCoefficientFraction(stokes_coeff_parallel, coeff.value()));
}

void StokesCoefficientFraction::calc_stokes_coeff() const
{
  stokes_coeff.value() = stokes_coeff_parallel;
  for(int i = 0; i < stokes_coeff.rows(); ++i) {
    stokes_coeff(i, 1) = stokes_coeff(i, 1) * (1 - 2 * coeff(i));
    stokes_coeff(i, 2) = stokes_coeff(i, 2) * (1 - 2 * coeff(i));
  }
}
