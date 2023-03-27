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

blitz::Array<double, 1> compute_raman_sioris(double solar_zenith, double viewing_zenith, double scattering_angle, double albedo, bool do_upwelling, const blitz::Array<double, 1> &temperature_layers, const blitz::Array<double, 1>& air_number_density, const SpectralDomain &grid, const SpectralDomain &grid_out, const blitz::Array<double, 1> &solar_irradiance, const blitz::Array<double, 2> &total_optical_depth);

/****************************************************************//**
 Implements adding the effect of Raman scattering using code
 originally by Christopher Sioris to compute the effect. There
 is a retrievable scale factor that accounts for the effect of
 multiple scattering.

 Note a few limitations here:

 1. We propagate jacobians through, but we don't properly account for
    the jacobian of the computed raman parameters. These depend on the
    optical depths and temperature, so they should depend on jacobians
    of these values. We simply ignore this. 
    (We do include d raman_scaling / d scale_factor, but not
     d raman_scaling / d temperature or d raman_scaling / d optical_depth).
 2. The code uses a simple surface model with a single albedo value.
    The general Ground models we have for the surface are more 
    complicated then this, and in general we don't even directly have
    an albedo. We currently only work with a SpurrBrdfType of
    LAMBERTIAN, we can examine that logic if we need to support other
    ground types.
 3. Because Raman scattering is looking for the "cross talk" between
    neighboring spectral points, it needs to have the optical depths
    and solar data for a wider range than we calculate the Raman
    data for. This given by raman_edge_wavenumber, and it hard coded
    in the Fortran code.
 4. Internally the Raman code works with a fixed space spectral grid
    that is exactly 1 wavenumber spacing. But we pass the solar and
    optical depth data at a coarser resolution that is then 
    interpolated to this 1 wavenumber spacing grid (and we then
    interpolate the resulting raman calculated grid to the spectral
    samplings of the Spectrum we are applying it to). There isn't
    any particularly obvious grid to use for the solar/optical depth,
    so for now we just take this as an input. The existing muses-py
    uses the same grid that the forward model is calculated on
    (since presumably the optical depths have already been calculated
    on this grid), but there is no reason to constrain this. So for
    now we expose this, and logic for how to determine this can be
    done in python (a reasonable choice might be the forward model
    grid plus the raman_edge_wavenumber padding).
*******************************************************************/
class RamanSiorisEffect : virtual public SpectrumEffectImpBase,
                          public Observer<Pressure> {

public:
  RamanSiorisEffect(const SpectralDomain&
                    Solar_and_odepth_spec_domain,
                    double scale_factor,
                    int channel_index, 
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

  virtual void notify_update(const Pressure& pressure) { compute_temp_layers(pressure); };
  virtual void notify_update(const StateVector& Sv)
  {
    // Call base notify to handle SV stuff
    SubStateVectorObserver::notify_update(Sv);
  };

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
    return res;
  } 

protected:

  void apply_raman_effect(Spectrum& Spec, const blitz::Array<double, 1>& temp_layers, const double albedo) const;

private:

  void compute_temp_layers(const Pressure& pressure);
  double evaluate_albedo(double wn, int cindex) const;

  SpectralDomain solar_and_odepth_spec_domain_;
  blitz::Array<double, 1> solar_and_odepth_wn_grid_;

  int channel_index_;
  bool do_upwelling_;

  double solar_zenith_;
  double obs_zenith_;
  double relative_azimuth_;
  double scattering_angle_;

  mutable blitz::Array<double, 1> temperature_layers_;

  boost::shared_ptr<AtmosphereStandard> atmosphere_;
  boost::shared_ptr<SolarModel> solar_model_;

  RamanSiorisEffect() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
};

}
FP_EXPORT_KEY(RamanSiorisEffect);
#endif
