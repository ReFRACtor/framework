#include "dispersion_polynomial.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(DispersionPolynomial, SampleGrid)
.def(luabind::constructor<const blitz::Array<double, 1>&, 
                          const blitz::Array<bool, 1>&,
                          const std::string&,
                          const blitz::Array<double, 1>&,
                          const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
/// Units passed seperately. Mainly used internally for cloning. Use
/// otherwise is considered deprecated.
//-----------------------------------------------------------------------

DispersionPolynomial::DispersionPolynomial(const blitz::Array<double, 1>& Coeff, 
                                           const blitz::Array<bool, 1>& Used_flag,
                                           const Unit& Coeff_unit,
                                           const blitz::Array<double, 1>& Var_values,
                                           const std::string& Band_name)
: SubStateVectorArray<SampleGrid>(Coeff, Used_flag),
  variable_values_(Var_values),
  coeff_unit(Coeff_unit),
  band_name_(Band_name),
  spectral_index(Var_values.rows())
{
  initialize();
}

//-----------------------------------------------------------------------
/// Constructor.
/// Pass units by string for use by Lua based unit testing. Consider this 
/// constructor deprecated.
//-----------------------------------------------------------------------

DispersionPolynomial::DispersionPolynomial(const blitz::Array<double, 1>& Coeff, 
                                           const blitz::Array<bool, 1>& Used_flag,
                                           const std::string& Coeff_unit_name,
                                           const blitz::Array<double, 1>& Var_values,
                                           const std::string& Band_name)
: SubStateVectorArray<SampleGrid>(Coeff, Used_flag),
  variable_values_(Var_values),
  coeff_unit(Coeff_unit_name),
  band_name_(Band_name),
  spectral_index(Var_values.rows())
{
  initialize();
}

//-----------------------------------------------------------------------
/// Preferred constructor.
//-----------------------------------------------------------------------

DispersionPolynomial::DispersionPolynomial(const ArrayWithUnit<double, 1>& Coeff, 
                                           const blitz::Array<bool, 1>& Used_flag,
                                           const blitz::Array<double, 1>& Var_values,
                                           const std::string& Band_name)
: SubStateVectorArray<SampleGrid>(Coeff.value, Used_flag),
  variable_values_(Var_values),
  coeff_unit(Coeff.units),
  band_name_(Band_name),
  spectral_index(Var_values.rows())
{
  initialize();
}

// Initialize class internals
void DispersionPolynomial::initialize() {
  // Initialize spectral_index to have the values 1....n
  firstIndex i1;
  spectral_index = i1 + 1;
}
 
// See base class for description.
std::string DispersionPolynomial::state_vector_name_i(int i) const
{
  std::string res = "Instrument Dispersion " + band_name_;
  if(i == 0)
    res += " Offset";
  else if(i == 1)
    res += " Scale";
  else
    res += " Parm " + boost::lexical_cast<std::string>(i + 1);
  return res;
}

// See base class for description.
SpectralDomain
DispersionPolynomial::pixel_grid() const
{
  Poly1d spectral_poly = Poly1d(coeff, false);
  ArrayAd<double, 1> index_array_ad(variable_values_, coeff.number_variable());
  index_array_ad.jacobian() = 0;
  SpectralDomain sample_grid = SpectralDomain(spectral_poly(index_array_ad), spectral_index, coeff_unit);
  return sample_grid;
}

boost::shared_ptr<SampleGrid> DispersionPolynomial::clone() const
{
  return boost::shared_ptr<SampleGrid>
    (new DispersionPolynomial(coeff.value(), used_flag, coeff_unit, variable_values_, band_name_));
}

void DispersionPolynomial::print(std::ostream& Os) const 
{
  Os << "DispersionPolynomial for band " << band_name_ << "\n"
     << "  Coeff:    (";
  for(int i = 0; i < coeff.rows() - 1; ++i)
    Os << coeff.value()(i) << ", ";
  Os << coeff.value()(coeff.rows() - 1) << ")\n"
     << "  Retrieve: (";
  for(int i = 0; i < used_flag.rows() - 1; ++i)
    Os << (used_flag(i) ? "true" : "false") << ", ";
  Os << (used_flag(used_flag.rows() - 1) ? "true" : "false") << ")";
}
