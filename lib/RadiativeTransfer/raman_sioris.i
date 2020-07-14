// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "raman_sioris.h"

// These are needed because of the inclusion of AtmosphereStandard
#include "instrument.h"
#include "spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
%}

%base_import(pressure)
%base_import(spectrum_effect_imp_base)

%import "spectral_domain.i"
%import "double_with_unit.i"
%import "atmosphere_standard.i"
%import "solar_model.i"

%fp_shared_ptr(FullPhysics::RamanSiorisEffect);
namespace FullPhysics {
  blitz::Array<double, 1> compute_raman_sioris(double solar_zenith,
       double viewing_zenith, double scattering_angle,
       double albedo, bool do_upwelling,
       const blitz::Array<double, 1> &temperature_layers,
       const blitz::Array<double, 1>& air_number_density,
       const SpectralDomain &grid,
       const blitz::Array<double, 1> &solar_irradiance,
       const blitz::Array<double, 2> &total_optical_depth);

class RamanSiorisEffect : public SpectrumEffectImpBase,
                          public Observer<Pressure> {

public:
    RamanSiorisEffect(double scale_factor, bool used_flag, 
                      int channel_index, 
                      const DoubleWithUnit& solar_zenith, 
                      const DoubleWithUnit& observation_zenith, 
                      const DoubleWithUnit& relative_azimuth,
                      const boost::shared_ptr<AtmosphereStandard>& atmosphere, 
                      const boost::shared_ptr<SolarModel>& solar_model,
                      double albedo,
                      double padding_fraction = 0.10,
                      bool do_upwelling = true,
                      double jac_perturbation = 0.001);

  virtual void apply_effect(Spectrum& Spec,
	      const ForwardModelSpectralGrid& Forward_model_grid) const;
  virtual void notify_update(const Pressure& pressure);
  virtual boost::shared_ptr<SpectrumEffect> clone() const;
  %pickle_serialization();
};
  
}

