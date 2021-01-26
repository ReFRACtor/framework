// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "ground.h"
#include "pressure.h"
%}

%base_import(observer)
%base_import(state_vector)
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::Ground)
namespace FullPhysics {
  class Ground;
}
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::Ground>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::Ground>)

namespace FullPhysics {
%template(ObservableGround) FullPhysics::Observable<FullPhysics::Ground>;
%template(ObserverGround) FullPhysics::Observer<FullPhysics::Ground>;

class Ground : public Observable<Ground> {
public:
  virtual ~Ground();
  virtual void add_observer(Observer<Ground>& Obs); 
  virtual void remove_observer(Observer<Ground>& Obs);
  std::string print_to_string() const;
  virtual ArrayAd<double, 1> surface_parameter
    (const double wn, const int spec_index) const = 0;
  virtual boost::shared_ptr<Ground> clone() const = 0;
  virtual void print(std::ostream& Os) const = 0;
  %pickle_serialization();
};
}

