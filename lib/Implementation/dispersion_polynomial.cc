#include "dispersion_polynomial.h"
#include "fp_serialize_support.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void DispersionPolynomial::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArraySampleGrid)
    & FP_NVP(coeff_unit) & FP_NVP_(band_name) & FP_NVP_(variable_values)
    & FP_NVP(spectral_index);
}

FP_IMPLEMENT(DispersionPolynomial);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
#include "state_mapping_at_indexes.h"

// Emulate old constructor that included flags by using a mapping class
// Return object at most generic level for Lua
boost::shared_ptr<SampleGrid> disp_poly_flagged_create(const blitz::Array<double, 1>& Coeff, 
                                                       const blitz::Array<bool, 1> Flag,
                                                       const std::string& Coeff_unit_name,
                                                       const blitz::Array<double, 1>& Var_values,
                                                       const std::string& Band_name)
{
    boost::shared_ptr<StateMapping> mapping =
        boost::make_shared<StateMappingAtIndexes>(Flag);

    boost::shared_ptr<DispersionPolynomial> disp_poly = 
        boost::make_shared<DispersionPolynomial>(Coeff, Coeff_unit_name, Var_values, Band_name, mapping);

    return disp_poly;
}

REGISTER_LUA_DERIVED_CLASS(DispersionPolynomial, SampleGrid)
.def(luabind::constructor<const blitz::Array<double, 1>&, 
                          const std::string&,
                          const blitz::Array<double, 1>&,
                          const std::string&>())
.scope
[
 luabind::def("create", &disp_poly_flagged_create)
]
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
/// Units passed seperately. Mainly used internally for cloning. Use
/// otherwise is considered deprecated.
//-----------------------------------------------------------------------

DispersionPolynomial::DispersionPolynomial(const blitz::Array<double, 1>& Coeff, 
                                           const Unit& Coeff_unit,
                                           const blitz::Array<double, 1>& Var_values,
                                           const std::string& Band_name,
                                           boost::shared_ptr<StateMapping> UNUSED(Mapping))
: coeff_unit(Coeff_unit),
  band_name_(Band_name),
  variable_values_(Var_values),
  spectral_index(Var_values.rows())
{
  SubStateVectorArray<SampleGrid>::init(Coeff);
  initialize();
}

//-----------------------------------------------------------------------
/// Constructor.
/// Pass units by string for use by Lua based unit testing. Consider this 
/// constructor deprecated.
//-----------------------------------------------------------------------

DispersionPolynomial::DispersionPolynomial(const blitz::Array<double, 1>& Coeff, 
                                           const std::string& Coeff_unit_name,
                                           const blitz::Array<double, 1>& Var_values,
                                           const std::string& Band_name,
                                           boost::shared_ptr<StateMapping> Mapping)

: coeff_unit(Coeff_unit_name),
  band_name_(Band_name),
  variable_values_(Var_values),
  spectral_index(Var_values.rows())
{
  SubStateVectorArray<SampleGrid>::init(Coeff, Mapping);
  initialize();
}

//-----------------------------------------------------------------------
/// Preferred constructor.
//-----------------------------------------------------------------------

DispersionPolynomial::DispersionPolynomial(const ArrayWithUnit<double, 1>& Coeff, 
                                           const blitz::Array<double, 1>& Var_values,
                                           const std::string& Band_name,
                                           boost::shared_ptr<StateMapping> Mapping)
: coeff_unit(Coeff.units),
  band_name_(Band_name),
  variable_values_(Var_values),
  spectral_index(Var_values.rows())
{
  SubStateVectorArray<SampleGrid>::init(Coeff.value, Mapping);
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
  ArrayAd<double, 1> poly_vals = mapping->mapped_state(coeff);
  Poly1d spectral_poly = Poly1d(poly_vals, false);
  ArrayAd<double, 1> index_array_ad(variable_values_, poly_vals.number_variable());
  index_array_ad.jacobian() = 0;
  SpectralDomain sample_grid = SpectralDomain(spectral_poly(index_array_ad), spectral_index, coeff_unit);
  return sample_grid;
}

boost::shared_ptr<SampleGrid> DispersionPolynomial::clone() const
{
  return boost::shared_ptr<SampleGrid>
    (new DispersionPolynomial(mapping->mapped_state(coeff).value(), coeff_unit, variable_values_, band_name_));
}

void DispersionPolynomial::print(std::ostream& Os) const 
{
  Os << "DispersionPolynomial for band " << band_name_ << "\n"
     << "  Coeff:    (";
  ArrayAd<double, 1> poly_vals = mapping->mapped_state(coeff);
  for(int i = 0; i < poly_vals.rows() - 1; ++i)
    Os << poly_vals.value()(i) << ", ";
  Os << poly_vals.value()(coeff.rows() - 1) << ")\n";
}
