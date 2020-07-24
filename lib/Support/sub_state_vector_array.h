#ifndef SUB_STATE_VECTOR_ARRAY_H
#define SUB_STATE_VECTOR_ARRAY_H
#include "sub_state_vector_observer.h"
#include "state_mapping_linear.h"
#include "pressure.h"
#include <blitz/array.h>
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>

namespace FullPhysics {
/****************************************************************//**
  It is common to have a class that is an Observable with a set of
  coefficients, a subset of which are updated by changes to the
  StateVector. This class captures this common behavior.

  For some classes, we only use the set of levels that lie above the
  lowest pressure level. For those classes, you can also pass in the
  Pressure class to use. If this doesn't apply to your particular
  class, then just leave the pressure out. This applies to the old
  fixed level classes (pre B2.10), this isn't currently used in the
  production code  anymore.
*******************************************************************/
template<class Base> class SubStateVectorArray:
    virtual public Base,
    virtual public SubStateVectorObserver {
public:
  // Note, we use to have constructors for SubStateVectorArray that
  // turn around and call init. Turns out this was causing a problem,
  // since SubStateVectorArray is generally inherited virtual (which
  // is required by boost serialization).
  //
  // This is because of a fairly subtle language rule. In C++ with
  // virtual base classes, the most derived class calls the
  // constructor. So if we have SubStateVectorArray<Foo> -> FooImpBase
  // -> FooExample then only FooExample calls SubStateVectorArray<Foo>
  // constructor. Most of the time this would means
  // SubStateVectorArray *default* constructor will likely be called
  // (see
  // https://stackoverflow.com/questions/9907722/why-is-default-constructor-called-in-virtual-inheritance).
  // This is almost certainly not what is intended, so we instead
  // don't depend on the constructor and instead call the init
  // function.
  //
  // We've removed all the nondefault constructors to keep from
  // accidentally calling them. Instead, use "init".
  //-----------------------------------------------------------------------
  /// Default constructor, should call init if we have Coeff to set.
  //-----------------------------------------------------------------------
  SubStateVectorArray()
  {
    // Default to no Coeff.
    blitz::Array<double, 1> no_coeff;
    blitz::Array<bool, 1> used_flag;
    init(no_coeff, used_flag);
  }

  void init
  (const blitz::Array<double, 1>& Coeff,
   const blitz::Array<bool, 1>& Used_flag,
   const boost::shared_ptr<Pressure>& Press = boost::shared_ptr<Pressure>(),
   bool Mark_according_to_press = true,
   int Pdep_start = 0,
   boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>())
  {
    mark_according_to_press = Mark_according_to_press;
    pdep_start = Pdep_start;
    if(Coeff.rows() != 0)
      coeff.reference(in_map->retrieval_init(Coeff.copy()));
    press = Press;
    if(Used_flag.rows() != 0)
      used_flag.reference(Used_flag.copy());
    mapping = in_map;
    
    if(coeff.rows() != used_flag.rows()) {
      throw Exception("Coeff and Used_flag need to be the same size");
    }
    
    state_vector_observer_initialize(count(Used_flag));
  }

  
  
  //-----------------------------------------------------------------------
  /// Special case when Coeff and Used_flag have exactly one row.
  //-----------------------------------------------------------------------
  void init(double Coeff, bool Used_flag,
    const boost::shared_ptr<Pressure>& Press = boost::shared_ptr<Pressure>(),
    bool Mark_according_to_press = true,
    int Pdep_start = 0,
    boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>())
  {
    blitz::Array<double, 1> coeff_t(1);
    blitz::Array<bool, 1> used_t(1);
    coeff_t(0) = Coeff;
    used_t(0) = Used_flag;
    init(coeff_t, used_t, Press, Mark_according_to_press, Pdep_start, in_map);
  }  
  virtual ~SubStateVectorArray() {}
  void mark_used_sub(blitz::Array<bool, 1>& Used) const
  {
    int si = 0;
    
    for(int i = 0; i < used_flag.rows(); ++i) {
      if(used_flag(i)) {
        if(!press || !mark_according_to_press ||
           i < press->number_level() + pdep_start) {
          Used(si) = true;
        }
        
        ++si;
      }
    }
  }


  //-----------------------------------------------------------------------
  /// Return a string to identify this part of the state, this name should be
  /// all lower case and seperate parts with a /. For example, an aerosol
  /// named strat would be named as:
  /// aerosol/strat.
  /// A gas named CO2 would be named like this:
  /// absorber/co2
  /// The name is intended to be used for looking up retrieval values 
  /// for a configuration system. Classes that have the same type of inputs
  /// should have the same name.
  //-----------------------------------------------------------------------
  
  virtual std::string sub_state_identifier() const
  {
    return "unknown/not_set";
  }
  
  //-----------------------------------------------------------------------
  /// Return state vector name for ith entry in coeff.
  //-----------------------------------------------------------------------
  
  virtual std::string state_vector_name_i(int i) const
  {
    return "Coeff" + boost::lexical_cast<std::string>(i + 1);
  }
  
  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) const
  {
    int si = 0;
        
    for(int i = 0; i < used_flag.rows(); ++i) {
      if(used_flag(i)) {
                Sv_name(si) = state_vector_name_i(i);
                ++si;
      }
    }
  }
  
  virtual void update_sub_state(const ArrayAd<double, 1>& Sv_sub, const blitz::Array<double, 2>& Cov)
  {
    if (Sv_sub.rows() > 0) {
      cov.reference(Cov.copy());
      int si = 0;
      coeff.resize_number_variable(Sv_sub.number_variable());
      
     for(int i = 0; i < coeff.rows(); ++i) {
        if(used_flag(i)) {
          coeff(i) = Sv_sub(si);
          ++si;
        }
     }
    } else {
      cov.resize(0,0);
      coeff.resize_number_variable(0);
    }
    
    update_sub_state_hook();
    Observable<Base>::notify_update_do(*this);
  }

  //-----------------------------------------------------------------------
  /// Hook for anything a derived class needs to do after coefficient is
  /// updated and before notify_update. Default is nothing.
  //-----------------------------------------------------------------------
  virtual void update_sub_state_hook()
  {
  }
  
  const ArrayAd<double, 1>& coefficient() const
  {
    return coeff;
  }
  
  const blitz::Array<bool, 1>& used_flag_value() const
  {
    return used_flag;
  }
  
  const blitz::Array<double, 2>& statevector_covariance() const
  {
    return cov;
  }
  
  const boost::shared_ptr<Pressure>& pressure() const
  {
    return press;
  }
protected:
  
  //-----------------------------------------------------------------------
  /// Coefficients.
  //-----------------------------------------------------------------------
  
  ArrayAd<double, 1> coeff;
  
  //-----------------------------------------------------------------------
  /// Pressure. This may be a null pointer, which just means this particular
  /// class doesn't store the Pressure object.
  //-----------------------------------------------------------------------
  
  boost::shared_ptr<Pressure> press;
  
  //-----------------------------------------------------------------------
  /// Flag indicating which of the coefficients gets updated by the
  /// StateVector.
  //-----------------------------------------------------------------------
  
  blitz::Array<bool, 1> used_flag;
  
  //-----------------------------------------------------------------------
  /// Last covariance matrix updated from the StateVector. If we haven't
  /// updated yet, this will be a 0x0 array.
  //-----------------------------------------------------------------------
  blitz::Array<double, 2> cov;
  
  //-----------------------------------------------------------------------
  /// Flag indicating if we only mark coefficients
  /// >= pdep_start + press->number_level() in mark_used_sub. This may seem
  /// a bit arcane, but this matches some of the classes that derive from
  /// this one (e.g., TemperatureFixedLevel).
  //-----------------------------------------------------------------------
  
  bool mark_according_to_press;
  
  //-----------------------------------------------------------------------
  /// Index of first coefficient that depends on the number of pressure
  /// levels. This is only used if mark_according_to_press to true, otherwise
  /// we don't do anything with this value.
  //-----------------------------------------------------------------------
  
  int pdep_start;
        
  //-----------------------------------------------------------------------
  /// StateMapping from internal coefficients to coefficients exposed to retrieval
  //-----------------------------------------------------------------------
  boost::shared_ptr<StateMapping> mapping;
        
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

#define SUB_STATE_VECTOR_ARRAY_SERIALIZE(Base, Type) \
template<> template<class Archive> \
void SubStateVectorArray<Base>::serialize(Archive & ar, \
                                          const unsigned int UNUSED(version)) \
{ \
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorObserver) \
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Base) \
    & FP_NVP(coeff) & FP_NVP(press) & FP_NVP(used_flag) & FP_NVP(cov) \
    & FP_NVP(mark_according_to_press) & FP_NVP(pdep_start) & FP_NVP(mapping); \
} \
\
FP_IMPLEMENT(Type);

#endif
