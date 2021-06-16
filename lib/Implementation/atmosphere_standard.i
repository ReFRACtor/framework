%include "fp_common.i"
%{
#include "atmosphere_standard.h"
#include "sub_state_vector_array.h"
#include "rayleigh.h"
#include "altitude.h"
%}
%fp_shared_ptr(FullPhysics::AtmosphereStandard);
%base_import(rt_atmosphere)
%base_import(aerosol)
%import "absorber.i"
%import "temperature.i"
%import "rayleigh.i"
%import "ground.i"
%import "surface_temperature.i"
%import "relative_humidity.i"
%import "altitude.i"
%import "constant.i"
namespace FullPhysics {
class AtmosphereStandard;
}

%fp_shared_ptr(FullPhysics::Observer<FullPhysics::AtmosphereStandard>);
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::AtmosphereStandard>);

namespace FullPhysics {
%template(ObserverAtmosphereStandard) FullPhysics::Observer<AtmosphereStandard>;

class AtmosphereStandard : public RtAtmosphere,
  public Observer<Aerosol>, public Observer<Pressure>,
  public CacheInvalidatedObserver {
public:
    // Supply all atmospheric consituent classes
    AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                       const boost::shared_ptr<Pressure>& pressurev,
                       const boost::shared_ptr<Temperature>& temperaturev,
                       const boost::shared_ptr<Rayleigh>& rayleighv,
                       const boost::shared_ptr<Aerosol>& aerosolv,
                       const boost::shared_ptr<RelativeHumidity>& rhv,
                       const boost::shared_ptr<Ground>& groundv,
                       const boost::shared_ptr<SurfaceTemperature>& surface_tempv,
                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                       const boost::shared_ptr<Constant>& C);

    // No surface temperature
    AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                       const boost::shared_ptr<Pressure>& pressurev,
                       const boost::shared_ptr<Temperature>& temperaturev,
                       const boost::shared_ptr<Rayleigh>& rayleighv,
                       const boost::shared_ptr<Aerosol>& aerosolv,
                       const boost::shared_ptr<RelativeHumidity>& rhv,
                       const boost::shared_ptr<Ground>& groundv,
                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                       const boost::shared_ptr<Constant>& C);

    // No ground, no surface temperature
    AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                       const boost::shared_ptr<Pressure>& pressurev,
                       const boost::shared_ptr<Temperature>& temperaturev,
                       const boost::shared_ptr<Rayleigh>& rayleighv,
                       const boost::shared_ptr<Aerosol>& aerosolv,
                       const boost::shared_ptr<RelativeHumidity>& rhv,
                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                       const boost::shared_ptr<Constant>& C);

    // No aerosol
    AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                       const boost::shared_ptr<Pressure>& pressurev,
                       const boost::shared_ptr<Temperature>& temperaturev,
                       const boost::shared_ptr<Rayleigh>& rayleighv,
                       const boost::shared_ptr<RelativeHumidity>& rhv,
                       const boost::shared_ptr<Ground>& groundv,
                       const boost::shared_ptr<SurfaceTemperature>& surface_tempv,
                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                       const boost::shared_ptr<Constant>& C);

    // No aerosol, no surface temperature
    AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                       const boost::shared_ptr<Pressure>& pressurev,
                       const boost::shared_ptr<Temperature>& temperaturev,
                       const boost::shared_ptr<Rayleigh>& rayleighv,
                       const boost::shared_ptr<RelativeHumidity>& rhv,
                       const boost::shared_ptr<Ground>& groundv,
                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                       const boost::shared_ptr<Constant>& C);

    // No ground, aerosol or surface temperature
    AtmosphereStandard(const boost::shared_ptr<Absorber>& absorberv,
                       const boost::shared_ptr<Pressure>& pressurev,
                       const boost::shared_ptr<Temperature>& temperaturev,
                       const boost::shared_ptr<Rayleigh>& rayleighv,
                       const boost::shared_ptr<RelativeHumidity>& Rh,
                       const std::vector<boost::shared_ptr<Altitude> >& altv,
                       const boost::shared_ptr<Constant>& C);

  virtual ArrayAdWithUnit<double, 1> altitude(int spec_index) const;
  %python_attribute_derived(number_spectrometer, int)
  %python_attribute_derived(number_layer, int)
  virtual AutoDerivative<double> 
  column_optical_depth(double wn, int spec_index, const std::string& Gas_name) const;
  virtual ArrayAd<double, 1> 
  optical_depth_wrt_rt(double wn, int spec_index) const;
  virtual ArrayAd<double, 1> 
  single_scattering_albedo_wrt_rt(double wn, int spec_index) const;
  virtual ArrayAd<double, 3>
  phase_function_moments_wrt_rt(double wn, int spec_index, int nummom = -1, 
                           int numscat = -1) const;
  virtual ArrayAd<double, 1> atmosphere_blackbody(double wn, int spec_index) const;
  virtual AutoDerivative<double> surface_blackbody(double wn, int spec_index) const;
  virtual boost::shared_ptr<OpticalProperties> optical_properties(double wn, int spec_index) const;
  virtual blitz::Array<double, 3> intermediate_jacobian(double wn, int spec_index) const;
  %python_attribute_derived(ground, boost::shared_ptr<Ground>)
  virtual void reset_timer();
  %python_attribute_derived(timer_info, std::string)
  %python_attribute2(pressure, pressure_ptr, boost::shared_ptr<Pressure>)
  %python_attribute2(absorber, absorber_ptr, boost::shared_ptr<Absorber>)
  %python_attribute2(aerosol, aerosol_ptr, boost::shared_ptr<Aerosol>)
  void set_aerosol(boost::shared_ptr<Aerosol>& new_aerosol, StateVector& Sv);
  %python_attribute2(temperature, temperature_ptr, boost::shared_ptr<Temperature>)
  %python_attribute2(relative_humidity, relative_humidity_ptr, boost::shared_ptr<RelativeHumidity>)
  %python_attribute2(constant, constant_ptr, boost::shared_ptr<Constant>)
  %python_attribute2(rayleigh, rayleigh_ptr, boost::shared_ptr<Rayleigh>)
  %python_attribute2(altitude_obj, altitude_ptr, std::vector<boost::shared_ptr<Altitude> >)
  boost::shared_ptr<AtmosphereStandard> clone() const;
  %python_attribute(rayleigh_only_atmosphere, bool)
};
}
