#ifndef ABSORBER_VMR_IMP_BASE_H
#define ABSORBER_VMR_IMP_BASE_H
#include "absorber_vmr.h"
#include "sub_state_vector_array.h"
#include "state_mapping_linear.h"
#include <boost/function.hpp>

namespace FullPhysics {
/****************************************************************//**
  As a design principle, we have each base class with the absolutely
  minimum interface needed for use from the rest of the system. This
  allows us to support any future code that supports this minimum 
  interface.
  
  However, almost always you will want to derive from this class 
  instead. See PressureImpBase for a more complete discussion of this.
*******************************************************************/
class AbsorberVmrImpBase: virtual public SubStateVectorArray<AbsorberVmr> {
public:
  virtual ~AbsorberVmrImpBase() {}
  virtual std::string gas_name() const {return gas_name_;}
  virtual AutoDerivative<double> 
  volume_mixing_ratio(const AutoDerivative<double>& P) const
  { fill_cache(); return vmr(P); }
  virtual boost::shared_ptr<AbsorberVmr> clone() const = 0;
  virtual void update_sub_state_hook() 
  { cache_stale = true; }
  virtual void print(std::ostream& Os) const { Os << "AbsorberVmrImpBase"; }

  virtual blitz::Array<bool, 1> state_used() const  {
    if (state_vector_start_index() == -1) {
      throw Exception("Must attach SubStateVectorObserver to state vector before checking used state.");
    }

    blitz::Array<bool, 1> res(sv_full.rows());
    res = false;
    if(state_vector_start_index() < res.rows())
      res(blitz::Range(state_vector_start_index(), 
           state_vector_start_index() + sub_vector_size() - 1)) = true;
    return res;
  }
protected:
//-----------------------------------------------------------------------
/// If this is true, the recalculate the vmr the next time we
/// need it.
//-----------------------------------------------------------------------
  mutable bool cache_stale;

//-----------------------------------------------------------------------
/// The cached volumn mixing ration. This should be filled in by
/// derived classes when calc_vmr() is called.
//-----------------------------------------------------------------------
  mutable boost::function<AutoDerivative<double>(AutoDerivative<double>)> vmr;

//-----------------------------------------------------------------------
/// Derived classes should provide a function to fill in vmr when this is 
/// called.
//-----------------------------------------------------------------------
  virtual void calc_vmr() const = 0;

//-----------------------------------------------------------------------
/// Initialize object.
//-----------------------------------------------------------------------

  void init
  (const std::string Gas_name,
   const blitz::Array<double, 1>& Coeff, 
   const boost::shared_ptr<Pressure>& Mapped_Press,
   boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>())
  {
    cache_stale = true;
    gas_name_ = Gas_name;
    mapped_pressure = Mapped_Press;
    SubStateVectorArray<AbsorberVmr>::init(Coeff, in_map);
  }

//-----------------------------------------------------------------------
/// Default constructor, derived class should call init if they use this
/// constructor.
//-----------------------------------------------------------------------

  AbsorberVmrImpBase() : cache_stale(true) { }

//-----------------------------------------------------------------------
/// Constructor that sets the coefficient() values.
//-----------------------------------------------------------------------
  AbsorberVmrImpBase(const std::string& Gas_name,
                     const blitz::Array<double, 1>& Coeff, 
                     const boost::shared_ptr<Pressure>& Mapped_Press,
                     boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>())
  {
    init(Gas_name, Coeff, Mapped_Press, in_map);
  }
protected:
  boost::shared_ptr<Pressure> mapped_pressure;
private:
  void fill_cache() const
  {
    if(cache_stale)
      calc_vmr();
    cache_stale = false;
  }
  std::string gas_name_;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
typedef SubStateVectorArray<AbsorberVmr> SubStateVectorArrayAbsorberVmr;
}

FP_EXPORT_KEY(AbsorberVmrImpBase);
FP_EXPORT_KEY(SubStateVectorArrayAbsorberVmr);
#endif
