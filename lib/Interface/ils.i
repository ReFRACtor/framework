%include <std_vector.i>
%include "fp_common.i"

%{
#include "ils.h"
%}

%base_import(state_vector_observer)

%import "array_ad.i"
%import "spectral_domain.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::Ils);

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Ils>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Ils>)

%template(ObservableIls) FullPhysics::Observable<FullPhysics::Ils>;
%template(ObserverIls) FullPhysics::Observer<FullPhysics::Ils>;

namespace FullPhysics {

%feature("director") Ils;

class Ils : public StateVectorObserver, public Observable<Ils> {
public:
  virtual ~Ils();
  virtual void add_observer(Observer<Ils>& Obs);
  virtual void remove_observer(Observer<Ils>& Obs);
  std::string print_to_string() const;
  std::string print_parent() const;
  virtual blitz::Array<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const blitz::Array<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const = 0;
  virtual ArrayAd<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const ArrayAd<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const = 0;
  virtual boost::shared_ptr<Ils> clone() const = 0;

  // When implementing directors can not use %python_attribute here or else there will
  // be missing symbol problems in the director
  virtual SpectralDomain pixel_grid() const = 0;
  virtual DoubleWithUnit high_res_extension() const = 0;
  virtual void high_res_extension(const DoubleWithUnit& extension) = 0;
  %pickle_serialization();
};
}

%template(vector_ils) std::vector<boost::shared_ptr<FullPhysics::Ils> >;

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(ils, Ils)

// List of things "import *" will include
%python_export("Ils", "ObservableIls", "ObserverIls");
