#ifndef PRESSURE_IMP_BASE_H
#define PRESSURE_IMP_BASE_H
#include "pressure.h"
#include "sub_state_vector_array.h"
#include "calculation_cache.h"

namespace FullPhysics {
class PressureImpBase;
  
/****************************************************************//**
  Cache used by PressureImpBase
*******************************************************************/
class PressureImpBaseCache : public CalculationCache<PressureImpBase> {
public:
  PressureImpBaseCache() {}
  virtual ~PressureImpBaseCache() {}
  virtual void fill_cache(const PressureImpBase& P);
  ArrayAdWithUnit<double, 1> pgrid;
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
  instead. This provides additional functionality that you will almost
  always want. We support have a subset set of the full StateVector
  provide coefficients for this class, as well as caching the 
  calculation of the pressure levels so we only do the calculation when
  something has changed (e.g., the StateVector).

  The more general Pressure interface can be used to support unusual
  cases that don't match this implementation, for example wrapping 
  an existing third party library that doesn't mesh with how
  PressureImpBase sets things up. But most of the time, you'll want
  to derive from this class.
*******************************************************************/
class PressureImpBase: virtual public SubStateVectorArray<Pressure> {
public:
  virtual ~PressureImpBase() {}
  virtual ArrayAdWithUnit<double, 1> pressure_grid(Pressure::PressureGridType Gtype = Pressure::INCREASING_PRESSURE) const;
  virtual boost::shared_ptr<Pressure> clone() const = 0;
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
  virtual std::string desc() const { return "PressureImpBase"; }

  virtual TypePreference type_preference() const {return type_preference_;}
  using SubStateVectorArray<Pressure>::update_sub_state;
  using SubStateVectorArray<Pressure>::state_vector_name_sub;
protected:
  mutable PressureImpBaseCache cache;
  Pressure::TypePreference type_preference_;

//-----------------------------------------------------------------------
/// Derived classes should provide a function to fill in pgrid when this is 
/// called.
//-----------------------------------------------------------------------
  virtual void calc_pressure_grid() const = 0;

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------

  PressureImpBase() : type_preference_(Pressure::PREFER_INCREASING_PRESSURE) { }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() values.
//-----------------------------------------------------------------------
  PressureImpBase(const blitz::Array<double, 1>& Coeff,
		  Pressure::TypePreference Tpref =
		  Pressure::PREFER_INCREASING_PRESSURE)
    : type_preference_(Tpref)
  {
    SubStateVectorArray<Pressure>::init(Coeff);
  }
private:
  friend PressureImpBaseCache;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
typedef SubStateVectorArray<Pressure> SubStateVectorArrayPressure;
}

FP_EXPORT_KEY(PressureImpBaseCache);
FP_EXPORT_KEY(PressureImpBase);
FP_CLASS_VERSION(PressureImpBase, 1);
FP_EXPORT_KEY(SubStateVectorArrayPressure);
#endif
