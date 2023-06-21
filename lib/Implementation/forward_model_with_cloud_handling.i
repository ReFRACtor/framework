// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"

%{
#include "forward_model_with_cloud_handling.h"
#include "sub_state_vector_array.h"
%}

%base_import(observer)
%base_import(state_vector_observer)
%base_import(sub_state_vector_array)
%base_import(forward_model)
%base_import(named_spectrum)
%import "generic_object_with_cloud_handling.i"

%fp_shared_ptr(FullPhysics::ForwardModelWithCloudHandling);
%fp_shared_ptr(FullPhysics::CloudFractionFromState);
%fp_shared_ptr(FullPhysics::CloudFraction);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::CloudFraction>);
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::CloudFraction>);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::CloudFraction>);

namespace FullPhysics {
class CloudFraction;  
%template(ObservableCloudFraction) Observable<CloudFraction>;
%template(ObserverCloudFraction) Observer<CloudFraction>;

class CloudFraction : virtual public StateVectorObserver, 
  public Observable<CloudFraction> { 
public: 
  CloudFraction(double Cfrac, 
		const boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>()); 
  %python_attribute(cloud_fraction, AutoDerivative<double>); 
  virtual boost::shared_ptr<CloudFraction> clone() const = 0; 
  virtual void add_observer(Observer<CloudFraction>& Obs); 
  virtual void remove_observer(Observer<CloudFraction>& Obs); 
  %python_attribute(subobject_list, std::vector<boost::shared_ptr<GenericObject> >); 
  std::string print_to_string() const; 
  %pickle_serialization(); 
}; 

%template(SubStateVectorArrayCloudFraction) FullPhysics::SubStateVectorArray<CloudFraction>; 

class CloudFractionFromState: 
    virtual public SubStateVectorArray<CloudFraction> { 
public: 
  CloudFractionFromState(double Cfrac, 
			 const boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>()); 
  virtual boost::shared_ptr<CloudFraction> clone() const; 
  %python_attribute_with_set(cloud_fraction, AutoDerivative<double>); 
}; 

class ForwardModelWithCloudHandling : public ForwardModel,
public Observable<boost::shared_ptr<FullPhysics::NamedSpectrum> > {
public: 
  ForwardModelWithCloudHandling 
  (const boost::shared_ptr<ForwardModel>& Fmodel, 
   boost::shared_ptr<CloudFraction>& Cfrac, 
   const std::vector<boost::shared_ptr<GenericObjectWithCloudHandling> >& 
   Cloud_handling_obj = std::vector<boost::shared_ptr<GenericObjectWithCloudHandling> >()); 
  const std::vector<boost::shared_ptr<GenericObjectWithCloudHandling> >& 
  cloud_handling_vector() const; 
  void add_cloud_handling_object 
  (const boost::shared_ptr<GenericObjectWithCloudHandling> & Obj); 
  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false) const; 
  void set_do_cloud(bool do_cloud) const; 
  virtual void add_observer(Observer<boost::shared_ptr<FullPhysics::NamedSpectrum> >& Obs); 
  virtual void remove_observer(Observer<boost::shared_ptr<FullPhysics::NamedSpectrum> >& Obs);
  void notify_spectrum_update(const Spectrum& updated_spec, const std::string& spec_name, int sensor_index) const;
  %python_attribute(underlying_forward_model, boost::shared_ptr<ForwardModel>); 
  %python_attribute(cloud_fraction, boost::shared_ptr<CloudFraction>); 
  %python_attribute(subobject_list, std::vector<boost::shared_ptr<GenericObject> >); 
}; 
}
