%include "fp_common.i"

%{
#include "scale_cloud.h"
%}

%base_import(state_vector_observer)
%import "auto_derivative.i"

%fp_shared_ptr(FullPhysics::ScaleCloud)

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::ScaleCloud>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::ScaleCloud>)

%template(ObservableScaleCloud) FullPhysics::Observable<FullPhysics::ScaleCloud>;
%template(ObserverScaleCloud) FullPhysics::Observer<FullPhysics::ScaleCloud>;

namespace FullPhysics {
class ScaleCloud : virtual public StateVectorObserver, public Observable<ScaleCloud> {
public:
  virtual ~ScaleCloud() {}
  virtual void add_observer(Observer<ScaleCloud>& Obs);
  virtual void remove_observer(Observer<ScaleCloud>& Obs);
  virtual AutoDerivative<double> scale_cloud(int sensor_index) const = 0;
  virtual boost::shared_ptr<ScaleCloud> clone() const = 0;
  virtual std::string desc() const;
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();
};
}

// List of things "import *" will include
%python_export("ScaleCloud");


