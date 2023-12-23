%include "fp_common.i"

%{
#include "array_ad.h"
#include "observation_sv_imp_base.h"
%}

%base_import(state_vector)
%base_import(sub_state_vector_array)
%base_import(observation)
%base_import(observation_sv)
%base_import(state_mapping)
%base_import(state_mapping_linear)
%import "state_mapping_linear.i"
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::ObservationSvImpBase);
%fp_shared_ptr(FullPhysics::SubStateVectorArray2<FullPhysics::ObservationSv, FullPhysics::Observation>);
namespace FullPhysics {
%template(SubStateVectorArrayObservationSv) 
    FullPhysics::SubStateVectorArray2<ObservationSv, Observation>;

// Allow these classes to be derived from in Python.
%feature("director") ObservationSvImpBase;

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
  
class ObservationSvImpBase: public SubStateVectorArray2<ObservationSv, Observation> {
public:
  void init
  (const blitz::Array<double, 1>& Coeff,
   boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
  void init(double Coeff,
            boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
  virtual std::string state_vector_name_i(int i) const;
  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) const;
  virtual void update_sub_state(const FullPhysics::ArrayAd<double, 1>& Sv_sub, const blitz::Array<double, 2>& Cov);
  virtual void update_sub_state_hook();
  %python_attribute(coefficient, FullPhysics::ArrayAd<double, 1>);
  %python_attribute(mapped_state, FullPhysics::ArrayAd<double, 1>);
  %python_attribute(sub_state_vector_values, FullPhysics::ArrayAd<double, 1>);
  %python_attribute(state_mapping, boost::shared_ptr<StateMapping>);
  %python_attribute(statevector_covariance, blitz::Array<double, 2>);
  %python_attribute_abstract(num_channels, int);
  virtual std::string desc() const;
  virtual SpectralDomain spectral_domain(int sensor_index) const = 0;
  boost::optional<blitz::Range> stacked_pixel_range(int sensor_index) const;
  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false)
    const = 0;
  virtual SpectralDomain spectral_domain_all() const;
  virtual Spectrum radiance_all(bool skip_jacobian = false) const;
  virtual std::string desc() const;
  %python_attribute(state_used, blitz::Array<bool, 1>)
  %sub_state_virtual_func(Observation);
  %pickle_serialization();
protected:
  FullPhysics::ArrayAd<double, 1> coeff;
  blitz::Array<double, 2> cov;
  boost::shared_ptr<StateMapping> mapping;  
  ObservationSvImpBase() {}
  ObservationSvImpBase(const blitz::Array<double, 1>& Coeff,
                boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(observation_sv_imp_base, ObservationSvImpBase)

// List of things "import *" will include
%python_export("ObservationSvImpBase", "SubStateVectorArrayObservationSv");
