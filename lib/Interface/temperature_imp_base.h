#ifndef TEMPERATURE_IMP_BASE_H
#define TEMPERATURE_IMP_BASE_H
#include "temperature.h"
#include "sub_state_vector_array.h"
#include "calculation_cache.h"
#include <boost/function.hpp>

namespace FullPhysics {
class TemperatureImpBase;
class TemperatureImpBaseCache : public CalculationCache<TemperatureImpBase> {
public:
  TemperatureImpBaseCache() {}
  virtual ~TemperatureImpBaseCache() {}
  virtual void fill_cache(const TemperatureImpBase& T);
  
//-----------------------------------------------------------------------
/// The cached temperature grid. This should be filled in by derived classes
/// when calc_temperature_grid() is called. This should map pressure
/// in Pascal to Temperature in Kelvin.
//-----------------------------------------------------------------------
  boost::function<AutoDerivative<double>(AutoDerivative<double>)> tgrid;
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
  
/****************************************************************//**
  As a design principle, we have each base class with the absolutely
  minimum interface needed for use from the rest of the system. This
  allows us to support any future code that supports this minimum 
  interface.
  
  However, almost always you will want to derive from this class 
  instead. See PressureImpBase for a more complete discussion of this.
*******************************************************************/
class TemperatureImpBase: virtual public SubStateVectorArray<Temperature> {
public:
  virtual ~TemperatureImpBase() {}
  virtual AutoDerivativeWithUnit<double> 
  temperature(const AutoDerivativeWithUnit<double>& Press) const
  { cache.fill_cache_if_needed(*this);
    return AutoDerivativeWithUnit<double>(cache.tgrid(Press.convert(units::Pa).value),
                                          units::K); 
  }
  virtual boost::shared_ptr<Temperature> clone() const = 0;
  
  virtual void update_sub_state_hook() 
  { cache.invalidate_cache(); }
  
//-----------------------------------------------------------------------
/// Print to stream. The default calls the function "desc" that returns
/// a string. This gives cleaner interface for deriving from this class
/// in python, but most C++ classes will want to override this function
/// rather than using desc.
//-----------------------------------------------------------------------
  virtual void print(std::ostream& Os) const { Os << desc(); }

//-----------------------------------------------------------------------
/// Description of object, to be printed to stream. This gives a cleaner
/// interface for deriving from python.
//-----------------------------------------------------------------------
  virtual std::string desc() const { return "TemperatureImpBase"; }
protected:
  friend TemperatureImpBaseCache;
  mutable TemperatureImpBaseCache cache;
  
//-----------------------------------------------------------------------
/// Derived classes should provide a function to fill in tgrid when this is 
/// called.
//-----------------------------------------------------------------------

  virtual void calc_temperature_grid() const = 0;

//-----------------------------------------------------------------------
/// Initialize object.
//-----------------------------------------------------------------------

  void init(const blitz::Array<double, 1>& Coeff, 
            const blitz::Array<bool, 1>& Used_flag,
            const boost::shared_ptr<Pressure>& Press,
            boost::shared_ptr<StateMapping> Map = boost::make_shared<StateMappingLinear>())

  { 
    SubStateVectorArray<Temperature>::init(Coeff, Used_flag, Map);
    press = Press;
  }

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------

  TemperatureImpBase() { }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() and used_flag() values.
//-----------------------------------------------------------------------
  TemperatureImpBase(const blitz::Array<double, 1>& Coeff, 
                     const blitz::Array<bool, 1>& Used_flag,
                     const boost::shared_ptr<Pressure>& Press,
                     boost::shared_ptr<StateMapping> Map = boost::make_shared<StateMappingLinear>())
  {
    init(Coeff, Used_flag, Press, Map);
  }

  boost::shared_ptr<Pressure> press;
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
typedef SubStateVectorArray<Temperature> SubStateVectorArrayTemperature;
}
FP_EXPORT_KEY(TemperatureImpBase);
FP_EXPORT_KEY(TemperatureImpBaseCache);
FP_EXPORT_KEY(SubStateVectorArrayTemperature);
#endif
