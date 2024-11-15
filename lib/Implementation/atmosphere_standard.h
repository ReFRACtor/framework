#ifndef ATMOSPHERE_STANDARD_H
#define ATMOSPHERE_STANDARD_H

#include "rt_atmosphere.h"

#include "constant.h"
#include "absorber.h"
#include "pressure.h"
#include "temperature.h"
#include "aerosol.h"
#include "ground.h"
#include "surface_temperature.h"
#include "relative_humidity.h"
#include "rayleigh.h"
#include "rayleigh_greek_moment.h"
#include "array_ad_cache.h"
#include "optical_properties_wrt_rt.h"
#include "altitude.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the atmosphere portion of the state, and uses
  this to set up the atmosphere and ground information needed to run
  the Radiative transfer code.

  This particular implementation forwards most of the work to other
  classes such as Absorber and Aerosol. This class then coordinates
  these other classes, and provides the calculations needed to set up
  the RT code.

  For some set ups, aerosol_ptr and/or ground_ptr may be null. For a
  Rayleigh only atmosphere, we don't have any aerosol to include. For
  up looking (e.g., TCCON FTS), there is no ground portion included in
  the radiative transfer.

  To speed up the calculation of the Jacobian in LIDORT, we make use
  of "intermediate" variables instead of directly using state vector
  variables. A description of this in more detail can be found in
  doc/LIDORT_Jacobian.pdf
*******************************************************************/
class AtmosphereStandard : public RtAtmosphere,
			   public Observer<Aerosol>,
			   public Observer<Pressure>,
			   public Observer<Absorber>,
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

  virtual ~AtmosphereStandard() {}

  virtual ArrayAdWithUnit<double, 1> altitude(int spec_index) const;

  virtual int number_spectrometer() const
  {
    return alt.size();
  }

  virtual int number_layer() const
  {
    if(nlay < 0) {
      nlay = pressure->number_layer();
    }

    return nlay;
  }

  virtual ArrayAd<double, 1> optical_depth_wrt_rt(double wn, int spec_index)
    const;

  virtual ArrayAd<double, 1> single_scattering_albedo_wrt_rt
  (double wn, int spec_index) const
  {
    fill_cache(wn, spec_index);
    return opt_prop->total_single_scattering_albedo();
  }

  virtual ArrayAd<double, 3> phase_function_moments_wrt_rt(double wn, int spec_index, int nummom = -1, int numscat = -1) const
  {
    fill_cache(wn, spec_index);
    return opt_prop->total_phase_function_moments(nummom, numscat);
  }

  virtual AutoDerivative<double> column_optical_depth
  (double wn, int spec_index, const std::string& Gas_name) const;

  virtual ArrayAd<double, 1> atmosphere_blackbody(double wn, int spec_index)
    const;

  virtual AutoDerivative<double> surface_blackbody(double wn, int spec_index)
    const;

  virtual const boost::shared_ptr<Ground> ground() const
  {
    return ground_ptr;
  }

  virtual boost::shared_ptr<OpticalProperties> optical_properties
  (double wn, int spec_index) const
  {
    fill_cache(wn, spec_index);

    return opt_prop;
  }

  // Use the state vector observer routines to update the length of the
  // state vector size used internally to allocated intermediate variables
  virtual void notify_add(StateVector& Sv)
  {
    sv_jac_size = (int) Sv.state_with_derivative().number_variable();

    // Can handle invalidation of across channel caches when SV changes
    can_cache_channel = true;
  }
  virtual void notify_remove(StateVector& Sv)
  {
    sv_jac_size = (int) Sv.state_with_derivative().number_variable();

    // Can no longer handle invalidation of across channel caches when
    // SV changes
    can_cache_channel = false;
  }
  virtual void notify_update(const StateVector& Sv)
  {
    notify_update_do(*this);
    sv_jac_size = (int) Sv.state_with_derivative().number_variable();

    // Invalidate caches since changes have been made to the atmosphere classes
    invalidate_local_cache();
  }

  virtual void print(std::ostream& Os) const;
  virtual void notify_update(const Aerosol& A);
  virtual void notify_update(const Absorber& UNUSED(A))
  {
    invalidate_local_cache();
    notify_update_do(*this);
  }
  virtual void notify_update(const Pressure& UNUSED(P))
  {
    if(nlay != pressure->number_layer()) {
      nlay = -1;
      invalidate_local_cache();
    }
  }
  virtual void invalidate_cache()
  {
    // Notify objects when *something* changed in AtmosphereStandard
    notify_update_do(*this);
  }
  virtual void reset_timer();
  virtual std::string timer_info() const;

  const boost::shared_ptr<Pressure>& pressure_ptr() const
  {
    return pressure;
  }
  const boost::shared_ptr<Absorber>& absorber_ptr() const
  {
    return absorber;
  }
  const boost::shared_ptr<Temperature>& temperature_ptr() const
  {
    return temperature;
  }
  const boost::shared_ptr<Aerosol>& aerosol_ptr() const
  {
    return aerosol;
  }
  void set_aerosol(boost::shared_ptr<Aerosol>& new_aerosol, StateVector& Sv);
  const boost::shared_ptr<RelativeHumidity>& relative_humidity_ptr() const
  {
    return rh;
  }
  const boost::shared_ptr<Constant>& constant_ptr() const
  {
    return constant;
  }
  const boost::shared_ptr<Rayleigh>& rayleigh_ptr() const
  {
    return rayleigh;
  }
  const std::vector<boost::shared_ptr<Altitude> >& altitude_ptr() const
  {
    return alt;
  }
  const boost::shared_ptr<Altitude>& altitude_ptr(int Spec_index) const
  {
    range_check(Spec_index, 0, number_spectrometer());
    return alt[Spec_index];
  }
  boost::shared_ptr<AtmosphereStandard> clone() const;

  //-----------------------------------------------------------------------
  /// Indicate we have rayleigh only atmosphere, i.e., we don't have any
  /// aerosol content.
  //-----------------------------------------------------------------------

  bool rayleigh_only_atmosphere() const
  {
    return !aerosol || aerosol->number_particle() == 0;
  }

  void attach_children_to_sv(StateVector& statev);
  virtual std::vector<boost::shared_ptr<GenericObject> >
  subobject_list() const
  { std::vector<boost::shared_ptr<GenericObject> > res;
    if(absorber)
      res.push_back(absorber);
    if(pressure)
      res.push_back(pressure);
    if(temperature)
      res.push_back(temperature);
    if(aerosol)
      res.push_back(aerosol);
    if(rh)
      res.push_back(rh);
    if(ground_ptr)
      res.push_back(ground_ptr);
    if(rayleigh)
      res.push_back(rayleigh);
    if(surface_temp)
      res.push_back(surface_temp);
    if(constant)
      res.push_back(constant);
    BOOST_FOREACH(auto i, alt)
      res.push_back(i);
    return res;
  }
private:

  boost::shared_ptr<Absorber> absorber;
  boost::shared_ptr<Pressure> pressure;
  boost::shared_ptr<Temperature> temperature;
  boost::shared_ptr<Aerosol> aerosol;
  boost::shared_ptr<RelativeHumidity> rh;
  boost::shared_ptr<Ground> ground_ptr;
  boost::shared_ptr<Rayleigh> rayleigh;
  boost::shared_ptr<SurfaceTemperature> surface_temp;
  boost::shared_ptr<Constant> constant;

  // The Altitude and Gravity constants depend on the specific
  // spectrometer we are using, because they see different ground
  // locations and this location enters into the hydrostatic equations
  // that determine the gravity and altitude constants.
  std::vector<boost::shared_ptr<Altitude> > alt;
  int sv_jac_size;

  // We cache the last calculation of these values
  // Keeping these around to reuse helps since
  // its almost certain each subsequent calculation
  // will have the same size, meaning the .resize operation
  // has no work to do and no new memory has to be allocated
  mutable int spec_index_tau_cache;
  mutable boost::shared_ptr<OpticalPropertiesWrtRt> opt_prop;
  mutable int nlay;

  // Items that might need to be cached for access
  // outside of say the radiative transfer loop, without
  // causing everything else to be recalculated
  boost::shared_ptr<ArrayAdCache<double, double, 1> > column_od_cache;
  boost::shared_ptr<ArrayAdCache<double, double, 1> > total_od_cache;
  mutable std::map<double, boost::shared_ptr<OpticalPropertiesWrtRt> > opt_prop_cache;
  mutable int n_cache_calls;
  mutable int n_cache_hits;
  mutable int n_cache_miss;

  // Only allow the column_od and total_od caches if the state vector
  // gets attached, otherwise the cache could become inconsistent with
  // the underlying sources of optical property information if values
  // are changed through an update through the state vector
  bool can_cache_channel;

  void initialize();
  bool fill_cache(double wn, int spec_index) const;
  void invalidate_local_cache() const;
  AtmosphereStandard();
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(AtmosphereStandard);
FP_CLASS_VERSION(AtmosphereStandard, 1);
#endif
