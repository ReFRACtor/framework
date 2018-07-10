#include "level_1b_sample_coefficient.h"
#include "fe_disable_exception.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA

// Convenience functions to give us spectral_coefficients in different forms for Lua
blitz::Array<double, 2> level_1b_s_coeffs(const Level1bSampleCoefficient& Lev1)
{
  int nparam = Lev1.spectral_coefficient(0).value.rows();
  blitz::Array<double, 2> res(Lev1.number_spectrometer(), nparam);
  for(int i = 0; i < res.rows(); ++i)
    res(i, blitz::Range::all()) = Lev1.spectral_coefficient(i).value;
  return res;
}

ArrayWithUnit<double, 2> level_1b_s_coeffs_with_unit(const Level1bSampleCoefficient& Lev1)
{
  int nparam = Lev1.spectral_coefficient(0).value.rows();
  Unit units = Lev1.spectral_coefficient(0).units;
  blitz::Array<double, 2> res(Lev1.number_spectrometer(), nparam);
  for(int i = 0; i < res.rows(); ++i)
    res(i, blitz::Range::all()) = Lev1.spectral_coefficient(i).value;
  return ArrayWithUnit<double,2>(res, units);
}


#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(Level1bSampleCoefficient, Level1b)
.def("spectral_coefficient", &level_1b_s_coeffs)
.def("spectral_coefficient_with_unit", &Level1bSampleCoefficient::spectral_coefficient)
.def("spectral_coefficient_with_unit", &level_1b_s_coeffs_with_unit)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes
typedef void(std::vector<boost::shared_ptr<Level1bSampleCoefficient> >::*pbt1)(
        const std::vector<boost::shared_ptr<Level1bSampleCoefficient> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<Level1bSampleCoefficient> >,
                        VectorLevel1bSampleCoefficient)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<Level1bSampleCoefficient> >::push_back))
REGISTER_LUA_END()

#endif

double Level1bSampleCoefficient::calculate_sample_value_from_coeffs(int Spec_index, int sample_index) const {
    ArrayWithUnit<double, 1> spectral_coefficients = this->spectral_coefficient(Spec_index);
    double sample_val = 0;
    for (int coeff_index = 0; coeff_index < spectral_coefficients.rows(); coeff_index++) {
        sample_val += (spectral_coefficients(coeff_index).value * std::pow(sample_index, coeff_index));
    }
    return sample_val;
}


SpectralDomain Level1bSampleCoefficient::sample_grid(int Spec_index) const {
    // Get number samples
    int num_samples = this->radiance(Spec_index).data().rows();
    Array<double, 1> grid_data = Array<double, 1>(num_samples);
    for (int sample_index = 0; sample_index < num_samples; sample_index++) {
        grid_data(sample_index) = this->calculate_sample_value_from_coeffs(Spec_index, sample_index);
    }
    ArrayWithUnit<double, 1> sample_grid_data = ArrayWithUnit<double, 1>(grid_data, this->spectral_coefficient(Spec_index).units);
    SpectralDomain sample_grid = SpectralDomain(sample_grid_data);
    return sample_grid;
}


