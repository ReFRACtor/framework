%include "common.i"

%{
#include "surface_temperature.h"
%}

%base_import(state_vector_observer)
%import "array_with_unit.i"
%import "auto_derivative_with_unit.i"

%fp_shared_ptr(FullPhysics::SurfaceTemperature)

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::SurfaceTemperature>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::SurfaceTemperature>)

%template(ObservableTemperature) FullPhysics::Observable<FullPhysics::SurfaceTemperature>;
%template(ObserverTemperature) FullPhysics::Observer<FullPhysics::SurfaceTemperature>;

namespace FullPhysics {
class SurfaceTemperature : virtual public StateVectorObserver, public Observable<SurfaceTemperature> {
public:
    virtual ~SurfaceTemperature() {}
    virtual void add_observer(Observer<SurfaceTemperature>& Obs);
    virtual void remove_observer(Observer<SurfaceTemperature>& Obs);
    virtual AutoDerivativeWithUnit<double> surface_temperature() const = 0;
    virtual boost::shared_ptr<SurfaceTemperature> clone() const = 0;
};
}


