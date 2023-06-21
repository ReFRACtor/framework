#ifndef FORWARD_MODEL_WITH_CLOUD_HANDLING_H
#define FORWARD_MODEL_WITH_CLOUD_HANDLING_H
#include "forward_model.h"
#include "sub_state_vector_array.h"
#include "generic_object_with_cloud_handling.h"
#include "named_spectrum.h"

namespace FullPhysics {
/****************************************************************//**
  Helper class to handle the cloud fraction. We could have just done
  this in the ForwardModelWithCloudHandling, but it seemed like
  a good idea to pull this out into a distinct class in case any
  other functionality needs to go here. 
*******************************************************************/

class CloudFraction : public Printable<CloudFraction>,
		      virtual public StateVectorObserver,
		      public Observable<CloudFraction>
{
public:
  CloudFraction() {}
  virtual ~CloudFraction() {}
  virtual void add_observer(Observer<CloudFraction>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<CloudFraction>& Obs) 
  { remove_observer_do(Obs, *this);}
  
//-----------------------------------------------------------------------
/// The cloud fraction.
//-----------------------------------------------------------------------

  virtual AutoDerivative<double> cloud_fraction() const = 0;

//-----------------------------------------------------------------------
/// Clone a CloudFraction object. Note that the cloned version
/// will *not* be attached to and StateVector or
/// Observer<CloudFraction>, although you can of course attach
/// them after receiving the cloned object. 
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" CloudFraction object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<CloudFraction> clone() const = 0;
  
//-----------------------------------------------------------------------
/// We have some fairly nested object hierarchies. It can be useful to
/// be able to search this for things (e.g., which Pressure object is
/// used by a ForwardModel?). This returns a list of subobjects
/// "owned" by this object.
//-----------------------------------------------------------------------

  virtual std::vector<boost::shared_ptr<GenericObject> >
  subobject_list() const
  { std::vector<boost::shared_ptr<GenericObject> > res;
    return res;
  }
  virtual void print(std::ostream& Os) const
  {
    Os << "Cloud Fraction: " << cloud_fraction();
  }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

typedef SubStateVectorArray<CloudFraction> SubStateVectorArrayCloudFraction;
  
class CloudFractionFromState:
    virtual public SubStateVectorArray<CloudFraction> {
public:
  CloudFractionFromState(double Cfrac,
		const boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>())
  {init(Cfrac, Mapping);}
  // TODO Fill this in
  virtual boost::shared_ptr<CloudFraction> clone() const
  { return boost::make_shared<CloudFractionFromState>(coefficient().value()(0),
						      mapping); }
  virtual ~CloudFractionFromState() {}
  virtual std::string sub_state_identifier() const 
  { return "cloud_fraction"; }
  virtual std::string state_vector_name_i(int UNUSED(i)) const
  { return "Cloud Fraction"; }
  virtual AutoDerivative<double> cloud_fraction() const { return coeff(0);}
  virtual void cloud_fraction(const AutoDerivative<double>& Cf) { coeff(0) = Cf;}
  virtual void print(std::ostream& Os) const
  {
    Os << "Cloud Fraction State Vector: " << cloud_fraction();
  }
private:
  CloudFractionFromState() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

/****************************************************************//**
  This is a forward model that handles clouds. We use a forward model
  to do both a clear and a fully cloudy retrieval, and then just 
  combine with a CloudFraction.

  The ForwardModel should contain objects like
  PressureWithCloudHandling and GroundWithCloudHandling. This class
  just takes a list of GenericObjectWithCloudHandling, and it toggles
  the "do_cloud" variable to get the clear and cloudy data.
*******************************************************************/

class ForwardModelWithCloudHandling : public ForwardModel,
		public Observable<boost::shared_ptr<NamedSpectrum> >{
public:
  ForwardModelWithCloudHandling
  (const boost::shared_ptr<ForwardModel>& Fmodel,
   const boost::shared_ptr<CloudFraction>& Cfrac,
   const std::vector<boost::shared_ptr<GenericObjectWithCloudHandling> >&
   Cloud_handling_vector = std::vector<boost::shared_ptr<GenericObjectWithCloudHandling> >());
  const std::vector<boost::shared_ptr<GenericObjectWithCloudHandling> >&
  cloud_handling_vector() const {return cloud_handling_vector_;}
  void add_cloud_handling_object
  (const boost::shared_ptr<GenericObjectWithCloudHandling> & Obj)
  { cloud_handling_vector_.push_back(Obj); }
  virtual ~ForwardModelWithCloudHandling() {}
  virtual void setup_grid() { fmodel_->setup_grid(); }
  virtual int num_channels() const { return fmodel_->num_channels(); }
  virtual SpectralDomain spectral_domain(int sensor_index) const
  { return fmodel_->spectral_domain(sensor_index); }
  virtual SpectralDomain::TypePreference spectral_domain_type_preference()
    const { return fmodel_->spectral_domain_type_preference(); }
  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false)
    const;
  virtual void print(std::ostream& Os) const;

//-----------------------------------------------------------------------
/// Change the clear/cloud state of the underlying forward model
/// (useful if you want to do a clear or cloudy only retrieval). Note
/// this is mostly for diagnostic sorts of things, you can do the
/// exact same thing by 
//-----------------------------------------------------------------------

  void set_do_cloud(bool do_cloud) const;
  
//-----------------------------------------------------------------------
/// The underlying forward model used for the clear and cloudy
/// retrievals.
//-----------------------------------------------------------------------

  const boost::shared_ptr<ForwardModel>& underlying_forward_model()
    const { return fmodel_; }

//-----------------------------------------------------------------------
/// The CloudFraction object.
//-----------------------------------------------------------------------

  const boost::shared_ptr<CloudFraction>& cloud_fraction() const
  { return cfrac_; }

  /// Required observable functions
  virtual void add_observer(Observer<boost::shared_ptr<NamedSpectrum> > & Obs)
  {
    add_observer_do(Obs);
  }

  virtual void remove_observer(Observer<boost::shared_ptr<NamedSpectrum> >& Obs)
  {
    remove_observer_do(Obs);
  }

  void notify_spectrum_update(const Spectrum& updated_spec, const std::string& spec_name, int sensor_index) const;

  //-----------------------------------------------------------------------
/// We have some fairly nested object hierarchies. It can be useful to
/// be able to search this for things (e.g., which Pressure object is
/// used by a ForwardModel?). This returns a list of subobjects
/// "owned" by this object.
//-----------------------------------------------------------------------

  virtual std::vector<boost::shared_ptr<GenericObject> >
  subobject_list() const
  { std::vector<boost::shared_ptr<GenericObject> > res;
    res.push_back(fmodel_);
    res.push_back(cfrac_);
    return res;
  }
private:
  boost::shared_ptr<ForwardModel> fmodel_;
  boost::shared_ptr<CloudFraction> cfrac_;
  std::vector<boost::shared_ptr<GenericObjectWithCloudHandling> >
  cloud_handling_vector_;
  ForwardModelWithCloudHandling() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(ForwardModelWithCloudHandling);
FP_EXPORT_KEY(CloudFraction);
FP_EXPORT_OBSERVER_KEY(CloudFraction);
FP_EXPORT_KEY(CloudFractionFromState);
FP_EXPORT_KEY(SubStateVectorArrayCloudFraction);
#endif
