// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "fp_common.i"

%{
#include "observation.h"
%}

%base_import(observer)
%base_import(observation)
%import "spectral_domain.i"
%import "spectrum.i"

%fp_shared_ptr(FullPhysics::Observation)
namespace FullPhysics {
  class Observation;
}
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Observation>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Observation>)

namespace FullPhysics {
%template(ObservableObservation) FullPhysics::Observable<FullPhysics::Observation>;
%template(ObserverObservation) FullPhysics::Observer<FullPhysics::Observation>;

class Observation : public Observable<Observation>
{
public:
  %python_attribute_abstract(num_channels, int)
  virtual void add_observer(Observer<Observation>& Obs); 
  virtual void remove_observer(Observer<Observation>& Obs);
  virtual SpectralDomain spectral_domain(int sensor_index,
				 bool include_bad_sample=false) const = 0;
  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false,
			    bool include_bad_sample=false) const = 0;
  virtual Spectrum radiance_all(bool skip_jacobian = false,
				bool include_bad_sample=false) const;
  %pickle_serialization();
};
}

%template(Vector_Observation) std::vector<boost::shared_ptr<FullPhysics::Observation> >;
