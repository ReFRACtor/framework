%include "fp_common.i"

%{
#include "absorber_vmr_imp_base.h"
%}
%base_import(sub_state_vector_array)
%base_import(state_vector)
%base_import(absorber_vmr)
%base_import(state_mapping)
%import "state_mapping_linear.i"

%fp_shared_ptr(FullPhysics::AbsorberVmrImpBase);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::AbsorberVmr>);
namespace FullPhysics {
%template(SubStateVectorArrayAbsorberVmr) 
     FullPhysics::SubStateVectorArray<AbsorberVmr>;

// Allow these classes to be derived from in Python.
%feature("director") AbsorberVmrImpBase;

// Note, at least for SWIG 2.0.4 a class that is derived from in python 
// needs to declare every virtual function that can be called on it, even 
// if all that happens is the base class to a director is called. This is 
// because this class is used to create the SwigDirector class, and this 
// class needs each of the member functions to direct things properly. It 
// is *not* necessary to add these function to the underlying
// C++, only that you declare them here.
//
// If you miss something, then you will get something like a recursion
// error in python when a virtual function is used that isn't explicitly
// listed here.
//
// This seems like a bug in 2.0.4, if SWIG needs all the member functions
// it should know to make them itself. So perhaps a future version of SWIG
// won't have this same constraint. But for now, this is required.

class AbsorberVmrImpBase: public SubStateVectorArray<AbsorberVmr> {
public:
  // From AbsorberVmrImpBase
  virtual ~AbsorberVmrImpBase();
  virtual std::string desc() const;
  virtual boost::shared_ptr<AbsorberVmr> clone() const = 0;
  %python_attribute(gas_name, std::string)
  virtual AutoDerivative<double> 
  volume_mixing_ratio(const AutoDerivative<double>& P) const;
  %python_attribute(state_used, blitz::Array<bool, 1>)

  %sub_state_virtual_func(AbsorberVmr);
  %pickle_serialization();
protected:
  mutable bool cache_stale;
  mutable boost::function<AutoDerivative<double>(AutoDerivative<double>)> vmr;
  virtual void calc_vmr() const = 0;
  AbsorberVmrImpBase();
  void init
  (const std::string Gas_name,
   const blitz::Array<double, 1>& Coeff, 
   const boost::shared_ptr<Pressure>& Mapped_Press);
  void init
  (const std::string Gas_name,
   const blitz::Array<double, 1>& Coeff, 
   const boost::shared_ptr<Pressure>& Mapped_Press,
   boost::shared_ptr<StateMapping> in_map);
  AbsorberVmrImpBase(const std::string& Gas_name,
                     const blitz::Array<double, 1>& Coeff,
                     const boost::shared_ptr<Pressure>& Press,
                     boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
  boost::shared_ptr<Pressure> mapped_pressure;
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(absorber_vmr_imp_base, AbsorberVmrImpBase)

// List of things "import *" will include
%python_export("AbsorberVmrImpBase");
