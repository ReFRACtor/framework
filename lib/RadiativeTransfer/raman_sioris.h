#ifndef RAMAN_SIORIS_H
#define RAMAN_SIORIS_H

#include <boost/shared_ptr.hpp>
#include <blitz/array.h>

#include "spectral_domain.h"
#include "spectrum_effect_imp_base.h"
#include "atmosphere_standard.h"
#include "solar_model.h"
#include "absorber.h"
#include "forward_model_spectral_grid.h"
#include "state_mapping_linear.h"

namespace FullPhysics {

blitz::Array<double, 1> compute_raman_sioris(double solar_zenith, double viewing_zenith, double scattering_angle, double albedo, bool do_upwelling, const blitz::Array<double, 1> &temperature_layers, const blitz::Array<double, 1>& air_number_density, const SpectralDomain &grid, const blitz::Array<double, 1> &solar_irradiance, const blitz::Array<double, 2> &total_optical_depth);

/****************************************************************//**
 Implements adding the effect of Raman scattering using code
 originally by Christopher Sioris to compute the effect. There
 is a retrievable scale factor that accounts for the effect of
 multiple scattering.
*******************************************************************/
class RamanSiorisEffect : virtual public SpectrumEffectImpBase,
                          public Observer<Pressure> {

public:
  RamanSiorisEffect(double scale_factor,
                    int channel_index, 
                    const DoubleWithUnit& solar_zenith, 
                    const DoubleWithUnit& observation_zenith, 
                    const DoubleWithUnit& relative_azimuth,
                    const boost::shared_ptr<AtmosphereStandard>& atmosphere, 
                    const boost::shared_ptr<SolarModel>& solar_model,
                    double albedo,
                    const boost::shared_ptr<StateMapping> mapping = boost::make_shared<StateMappingLinear>(),
                    double padding_fraction = 0.10,
                    bool do_upwelling = true,
                    double jac_perturbation = 0.001);

  virtual void apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& Forward_model_grid) const;

  virtual void notify_update(const Pressure& pressure) { compute_temp_layers(pressure); };

  virtual boost::shared_ptr<SpectrumEffect> clone() const;

  virtual std::string sub_state_identifier() const 
  { return "raman_sioris_" + boost::lexical_cast<std::string>(channel_index_ + 1); }

  virtual std::string state_vector_name_i(int UNUSED(i)) const
  { return "Raman Sioris Scale Factor, Channel #" + boost::lexical_cast<std::string>(channel_index_ + 1); }

  virtual void print(std::ostream& Os) const;

  virtual std::string name() const { return "raman_sioris"; }
  virtual std::vector<boost::shared_ptr<GenericObject> >
  subobject_list() const
  { std::vector<boost::shared_ptr<GenericObject> > res;
    res.push_back(atmosphere_);
    res.push_back(solar_model_);
    res.push_back(absorber_);
    return res;
  }

private:

  void compute_temp_layers(const Pressure& pressure);
  double evaluate_albedo() const;

  int channel_index_;
  double albedo_;
  double padding_fraction_;
  bool do_upwelling_;
  double jac_perturbation_;

  double solar_zenith_;
  double obs_zenith_;
  double relative_azimuth_;
  double scattering_angle_;

  blitz::Array<double, 1> temperature_layers_;

  boost::shared_ptr<AtmosphereStandard> atmosphere_;
  boost::shared_ptr<SolarModel> solar_model_;

  boost::shared_ptr<Absorber> absorber_;
  RamanSiorisEffect() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}
FP_EXPORT_KEY(RamanSiorisEffect);
#endif
