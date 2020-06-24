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

  The x values for the evaluated polynomical ranges from 0 to the
  number of samples. The offset constructor argument is a number
  added to each x value. Hence, x becomes x + offset.
*******************************************************************/
class Level1bSampleCoefficient: public Level1b {

public:

  virtual ~Level1bSampleCoefficient() = default;

//-----------------------------------------------------------------------
/// Number of samples of data corresponding to the sample_grid size for
//  a given instrument channel
//-----------------------------------------------------------------------

  virtual int number_sample(int channel_index) const = 0;

//-----------------------------------------------------------------------
/// Returns coefficients for an equation describing the special domain
/// used to translate radiance value indexes to their corresponding 
/// spectral grid. (ie wavenumber, wavelength, etc)
/// The meaning of these coefficients will be specific to the instrument
/// that measured the data.
//-----------------------------------------------------------------------

  virtual ArrayWithUnit<double, 1> spectral_coefficient(int channel_index) const = 0;

//-----------------------------------------------------------------------
/// Return the spectral coefficient variable values that should be 
/// evaluated at each sample index. In a simple case this might just be 
/// a zero based or one based indexing of the samples.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> spectral_variable(int channel_index) const = 0;

//-----------------------------------------------------------------------
/// Returns the sample grid (ie wavenumber, wavelength, etc) for the
/// corresponding radiance values.
/// The meaning of these coefficients will be specific to the instrument
/// that measured the data.
//-----------------------------------------------------------------------

  virtual SpectralDomain sample_grid(int channel_index) const;

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const {Os << "Level1bSampleCoefficient";}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
} // End of FullPhysics namespace

FP_EXPORT_KEY(Level1bSampleCoefficient);

#endif
