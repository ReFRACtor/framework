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

%feature("director") RamanSiorisEffect;

%fp_shared_ptr(FullPhysics::RamanSiorisEffect);
namespace FullPhysics {
  blitz::Array<double, 1> compute_raman_sioris(double solar_zenith,
       double viewing_zenith, double scattering_angle, double albedo,
       bool do_upwelling,
       const blitz::Array<double, 1> &temperature_layers,
       const blitz::Array<double, 1>& air_number_density,
       const SpectralDomain &grid,
       const SpectralDomain &grid_out,
       const blitz::Array<double, 1> &solar_irradiance,
       const blitz::Array<double, 2> &total_optical_depth);

class RamanSiorisEffect : public SpectrumEffectImpBase,
                          public Observer<Pressure> {

public:
  RamanSiorisEffect(const SpectralDomain& Solar_and_odepth_spec_domain,
                    double scale_factor,
                    int sensor_index, 
                    const DoubleWithUnit& solar_zenith, 
                    const DoubleWithUnit& observation_zenith, 
                    const DoubleWithUnit& relative_azimuth,
                    const boost::shared_ptr<AtmosphereStandard>& atmosphere, 
                    const boost::shared_ptr<SolarModel>& solar_model,
                    const boost::shared_ptr<StateMapping> mapping = boost::make_shared<StateMappingLinear>(),
                    bool do_upwelling = true);

  /// The "edge" we need to the desired range of the Raman calculation
  static const double raman_edge_wavenumber;
  
  virtual void apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& Forward_model_grid) const;

  virtual void notify_update(const Pressure& pressure);
  virtual void notify_update(const StateVector& Sv);

  virtual boost::shared_ptr<SpectrumEffect> clone() const;

  virtual std::string sub_state_identifier() const;
  virtual std::string state_vector_name_i(int i) const;

  // Indicate this is implemented by this class since defined in SpectrumEffect
  %python_attribute(name, std::string);
  %python_attribute(subobject_list, std::vector<boost::shared_ptr<GenericObject> >);

  virtual std::string desc() const;
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();

protected:
  RamanSiorisEffect();
  void apply_raman_effect(Spectrum& Spec, const blitz::Array<double, 1>& temp_layers, const double albedo) const;

};

const double RamanSiorisEffect::raman_edge_wavenumber = 218;


}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(raman_sioris, RamanSiorisEffect)

// List of things "import *" will include
%python_export("RamanSiorisEffect", "compute_raman_sioris");
