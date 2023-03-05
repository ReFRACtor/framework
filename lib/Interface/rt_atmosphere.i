// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "rt_atmosphere.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}

%base_import(observer)
%base_import(state_vector_observer)
%import "array_ad_with_unit.i"
%import "auto_derivative.i"
%import "ground.i"
%import "optical_properties.i"

%fp_shared_ptr(FullPhysics::RtAtmosphere);

namespace FullPhysics {
  class RtAtmosphere;
}

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::RtAtmosphere>);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::RtAtmosphere>);

namespace FullPhysics {
%template(ObservableRtAtmosphere) Observable<RtAtmosphere>;
%template(ObserverRtAtmosphere) Observer<RtAtmosphere>;

class RtAtmosphere : virtual public StateVectorObserver, 
                 public Observable<RtAtmosphere> {
public:
  virtual ~RtAtmosphere();
  virtual void add_observer(Observer<RtAtmosphere>& Obs);
  virtual void remove_observer(Observer<RtAtmosphere>& Obs);
  %python_attribute(timer_info, std::string)
  %python_attribute_abstract(number_layer, int)
  %python_attribute_abstract(number_spectrometer, int)
  virtual ArrayAdWithUnit<double, 1> altitude(int spec_index) const = 0;
  virtual AutoDerivative<double> 
  column_optical_depth(double wn, int spec_index, const std::string& Gas_name) const = 0;
  virtual ArrayAd<double, 1> optical_depth_wrt_rt(double wn, int spec_index) const = 0;
  virtual ArrayAd<double, 1> single_scattering_albedo_wrt_rt(double wn, int spec_index) const = 0;
  virtual ArrayAd<double, 3> phase_function_moments_wrt_rt(double wn, int spec_index, int nummom = -1, int numscat = -1) const = 0;
  virtual ArrayAd<double, 1> optical_depth_wrt_state_vector(double wn, int spec_index) const;
  virtual ArrayAd<double, 1> single_scattering_albedo_wrt_state_vector(double wn, int spec_index) const;
  virtual ArrayAd<double, 3> phase_function_moments_wrt_state_vector(double wn, int spec_index, int nummom = -1, int numscat = -1) const;
  virtual ArrayAd<double, 1> atmosphere_blackbody(double wn, int spec_index) const = 0;
  virtual AutoDerivative<double> surface_blackbody(double wn, int spec_index) const = 0;
  virtual boost::shared_ptr<OpticalProperties> optical_properties(double wn, int spec_index) const = 0;
  virtual blitz::Array<double, 3> intermediate_jacobian(double wn, int spec_index) const;
  %python_attribute_abstract(ground, boost::shared_ptr<Ground>)
  virtual void reset_timer();
  std::string print_to_string() const;
  %python_attribute(subobject_list, std::vector<boost::shared_ptr<GenericObject> >);
  %pickle_serialization()
};
}

