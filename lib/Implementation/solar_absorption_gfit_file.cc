#include "solar_absorption_gfit_file.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "linear_interpolate.h"

using namespace FullPhysics;
using namespace blitz;
#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SolarAbsorptionGfitFile::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SolarAbsorptionSpectrum)
    & FP_NVP_(line_list_file) & FP_NVP_(fraction_solar_diameter);
}

FP_IMPLEMENT(SolarAbsorptionGfitFile);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_DERIVED_CLASS(SolarAbsorptionGfitFile, SolarAbsorptionSpectrum)
.def(luabind::constructor<const std::string&, double>())
REGISTER_LUA_END()
#endif

extern "C" {
  void solar_pts(const int *lunr, const char* filename, const int* filename_len, const double *fzero, const double *grid, const double *frac, double *spts, const int *ncp);
}

//-----------------------------------------------------------------------
/// Read the given line list file, and use for calculating the solar
/// absorption spectrum.
///
/// \param Line_list_file Line list file
/// \param Fraction_solar_diameter Fraction of Solar diameter viewed.
//-----------------------------------------------------------------------

SolarAbsorptionGfitFile::SolarAbsorptionGfitFile(const std::string& Line_list_file,
                                                 double Fraction_solar_diameter)
: line_list_file_(Line_list_file), fraction_solar_diameter_(Fraction_solar_diameter)
{
}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

void SolarAbsorptionGfitFile::print(std::ostream& Os) const
{
  Os << "Solar Absorption GFIT File:\n"
     << "  Fraction solar diameter: " << fraction_solar_diameter_;
}

// See base class for description.
Spectrum SolarAbsorptionGfitFile::solar_absorption_spectrum( const SpectralDomain& spec_domain) const
{
  // Set up inputs to match our grid
  // Convert to wavenumber regardless of the input grid units,
  // Find the mean spacing which might vary if converted from wavenumbers
  // As long as we abs the difference it should not matter the order
  // Grid used inside of routine is computed as:
  //  V(i) = fzero + i * grid       i = 1,NCP
  Array<double, 1> wavenumbers = spec_domain.wavenumber();

  Array<double, 1> wn_1 = wavenumbers(Range(1, wavenumbers.rows()-1));
  Array<double, 1> wn_0 = wavenumbers(Range(0, wavenumbers.rows()-2));

  double grid = mean(abs(wn_1 - wn_0));
  double fzero = min(wavenumbers);
  int ncp = spec_domain.data().rows();

  int lunr = 99;
  int file_len = line_list_file_.length();

  Array<double, 1> spts_calc(ncp); // output array
  solar_pts(&lunr, line_list_file_.c_str(), &file_len, &fzero, &grid, &fraction_solar_diameter_, spts_calc.dataFirst(), &ncp);

  // Interpolate back to input grid
  Array<double, 1> spts_grid(ncp);
  for (int wn_idx = 0; wn_idx < ncp; wn_idx++) {
      spts_grid(wn_idx) = fzero + wn_idx * grid;
  }

  LinearInterpolate<double, double> spts_interp
    (spts_grid.dataFirst(), spts_grid.dataFirst() + spts_grid.rows(), spts_calc.dataFirst());

  Array<double, 1> spts_out(spec_domain.data().rows());
  for (int wn_idx = 0; wn_idx < spec_domain.data().rows(); wn_idx++) {
    spts_out(wn_idx) = spts_interp(wavenumbers(wn_idx));
  }

  return Spectrum(spec_domain, SpectralRange(spts_out, units::dimensionless));
}
