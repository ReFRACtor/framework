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
REGISTER_LUA_CLASS(Level1bSampleCoefficient)
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


/* TODO: Implement return of SpectralDomain from spectral_coefficient */
SpectralDomain Level1bSampleCoefficient::sample_spectral_domain(int Spec_index) const {
    // Get number samples
    SpectralRange spec_radiance = this->radiance(Spec_index);
    Array<double, 1> rad_data = spec_radiance.data();
    int number_samples = rad_data.extent(0);
    cout << "Hi! Got number_samples" << number_samples << endl;
    // Allocate ArraywithUnit with that number of samples
    // For each sample, calculate wavelength from polynomial
    // Create spectralDomain from data
    ArrayWithUnit<double, 1> data;
    SpectralDomain _temp = SpectralDomain(data);
    return _temp;
}

#endif
