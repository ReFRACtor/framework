%include "fp_common.i"

%{
#include "observation_sv_imp_base.h"
%}

%base_import(state_vector)
%base_import(sub_state_vector_array)
%base_import(observation)
%base_import(observation_sv)
%base_import(state_mapping)
%base_import(state_mapping_linear)
%import "state_mapping_linear.i"

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
  %python_attribute_abstract(num_channels, int);
  virtual SpectralDomain spectral_domain(int sensor_index, bool include_bad_sample=false) const = 0;
  boost::optional<blitz::Range> stacked_pixel_range(int sensor_index) const;
  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false, bool include_bad_sample=false)
    const = 0;
  virtual Spectrum radiance_all(bool skip_jacobian = false, bool include_bad_sample=false) const;
  virtual std::string desc() const;
  %python_attribute(state_used, blitz::Array<bool, 1>)
  %sub_state_virtual_func(Observation);
  %pickle_serialization();
protected:
  ObservationSvImpBase(const blitz::Array<double, 1>& Coeff,
                boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
};
}

