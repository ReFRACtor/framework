#ifndef LEVEL_1B_SAMPLE_COEFFICIENT_H
#define LEVEL_1B_SAMPLE_COEFFICIENT_H
#include "level_1b.h"
#include "spectral_domain.h"
#include "array_ad.h"
#include "polynomial_eval.h"
#include <blitz/array.h>
#include <stdint.h>
#include <cmath>

namespace FullPhysics {
/****************************************************************//**
  This is used to wrap the nominal Level 1B file reader.
  It creates the required SpectralDomain from the given
  coefficients.
  It calculautes the wavelength/wavenumber as follows:
  f(x) = x*(coeff[0]^0) + x*(coeff[1]^1) ... + x*(coeff[n]^n)
  Where n is spectral_coefficient length.
*******************************************************************/
class Level1bSampleCoefficient: public Level1b {
public:
  virtual ~Level1bSampleCoefficient() { };

//-----------------------------------------------------------------------
/// Returns coefficients for an equation describing the special domain
/// used to translate radiance value indexes to their corresponding 
/// spectral grid. (ie wavenumber, wavelength, etc)
/// The meaning of these coefficients will be specific to the instrument
/// that measured the data.
//-----------------------------------------------------------------------

  virtual ArrayWithUnit<double, 1> spectral_coefficient(int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// Returns the sample grid (ie wavenumber, wavelength, etc) for the
/// corresponding radiance values.
/// The meaning of these coefficients will be specific to the instrument
/// that measured the data.
//-----------------------------------------------------------------------

  virtual SpectralDomain sample_grid(int Spec_index) const;

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const {Os << "Level1bSampleCoefficient";}

private:
  virtual double calculate_sample_value_from_coeffs(int Spec_index, int sample_idx) const;
};
} // End of FullPhysics namespace
#endif
