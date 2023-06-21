// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "observation_sv.h"
%}

%base_import(observer)
%base_import(observation)
%base_import(state_vector)
%base_import(state_vector_observer)
%import "spectrum.i"

%fp_shared_ptr(FullPhysics::ObservationSv)
namespace FullPhysics {
  class ObservationSv;
}
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::ObservationSv>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::ObservationSv>)

namespace FullPhysics {
%template(ObservableObservationSv) FullPhysics::Observable<FullPhysics::ObservationSv>;
%template(ObserverObservationSv) FullPhysics::Observer<FullPhysics::ObservationSv>;
  
class ObservationSv : public Observable<ObservationSv>, public Observation,
		      public StateVectorObserver {
public:
  virtual ~ObservationSv();
  virtual void add_observer(Observer<ObservationSv>& Obs); 
  virtual void remove_observer(Observer<ObservationSv>& Obs);
  virtual boost::shared_ptr<ObservationSv> clone() const;
  %pickle_serialization();
};
}

