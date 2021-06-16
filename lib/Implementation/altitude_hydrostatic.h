#ifndef ALTITUDE_HYDROSTATIC_H
#define ALTITUDE_HYDROSTATIC_H
#include "altitude.h"
#include "pressure.h"
#include "temperature.h"
#include "double_with_unit.h"
#include "linear_interpolate.h"
#include "calculation_cache.h"

namespace FullPhysics {
class AltitudeHydrostatic;
class AltitudeHydrostaticCache : public CalculationCache<AltitudeHydrostatic>
{
public:
  AltitudeHydrostaticCache() {}
  virtual ~AltitudeHydrostaticCache() {}
  virtual void fill_cache(const AltitudeHydrostatic& A);
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >
    lin_type;
  boost::shared_ptr<lin_type> alt;
  boost::shared_ptr<lin_type> grav;
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
  
/****************************************************************//**
  This class handles the calculation of the altitude an gravity
  constants, automatically updating with the surface pressure or
  temperature profile is updated.

  We do this by solving the hydrostatic equations.

  \todo reference for this? (ATB?)
*******************************************************************/
class AltitudeHydrostatic : public Altitude
{
public:
  AltitudeHydrostatic(const boost::shared_ptr<Pressure>& P,
                      const boost::shared_ptr<Temperature>& T,
                      const DoubleWithUnit& Latitude, 
                      const DoubleWithUnit& Surface_height,
                      const int Num_sublayer = 10);
  virtual ~AltitudeHydrostatic() {};
  virtual AutoDerivativeWithUnit<double> altitude(const AutoDerivativeWithUnit<double>& P) const
  {
    cache->fill_cache_if_needed(*this);
    AutoDerivativeWithUnit<double> p_pas = P.convert(units::Pa);
    return AutoDerivativeWithUnit<double>((*cache->alt)(p_pas.value), units::km); 
  }

  virtual AutoDerivativeWithUnit<double> gravity(const AutoDerivativeWithUnit<double>& P) const
  {
    cache->fill_cache_if_needed(*this);
    AutoDerivativeWithUnit<double> p_pas = P.convert(units::Pa);
    return AutoDerivativeWithUnit<double>((*cache->grav)(p_pas.value), "m/s^2"); 
  }

  virtual void print(std::ostream& Os) const { Os << "AltitudeHydrostatic"; }
  virtual boost::shared_ptr<Altitude> clone() const;
private:
  DoubleWithUnit latitude, surface_height;
  friend AltitudeHydrostaticCache;
  mutable boost::shared_ptr<AltitudeHydrostaticCache> cache;
  boost::shared_ptr<Pressure> p;
  boost::shared_ptr<Temperature> t;
  int num_sublayer;
  AutoDerivative<double> 
  gravity_calc(double gdlat, AutoDerivative<double> altit) const;
  void altitude_calc
  (const ArrayAdWithUnit<double, 1>& press_grid,  
   const ArrayAdWithUnit<double, 1>& temp_grid) const;
  AltitudeHydrostatic() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(AltitudeHydrostaticCache);
FP_EXPORT_KEY(AltitudeHydrostatic);
#endif
