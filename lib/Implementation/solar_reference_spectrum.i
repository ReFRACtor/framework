// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"

%{
#include "solar_reference_spectrum.h"
#include "sub_state_vector_array.h"
#include "instrument.h"
#include "spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
%}

%base_import(solar_model)
%import "solar_doppler_shift.i"

%fp_shared_ptr(FullPhysics::SolarReferenceSpectrum);

namespace FullPhysics {
class SolarReferenceSpectrum : public SolarModel {
public:
  SolarReferenceSpectrum(const boost::shared_ptr<Spectrum>& reference_spectrum,
			 const boost::shared_ptr<SolarDopplerShift>& doppler_shift = NULL);

  virtual ~SolarReferenceSpectrum() = default;

  virtual boost::shared_ptr<SpectrumEffect> clone() const;

  virtual const boost::shared_ptr<SolarDopplerShift>& doppler_shift() const;

  virtual void print(std::ostream& Os) const;
  virtual Spectrum solar_spectrum(const SpectralDomain& spec_domain) const;
  %pickle_serialization();
};
}
