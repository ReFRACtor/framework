// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "altitude_hydrostatic.h"
%}
%base_import(altitude)
%import "double_with_unit.i"
%import "auto_derivative_with_unit.i"
%import "pressure.i"
%import "temperature.i"
%fp_shared_ptr(FullPhysics::AltitudeHydrostatic)
namespace FullPhysics {
class AltitudeHydrostatic : public Altitude {
public:
  AltitudeHydrostatic(const boost::shared_ptr<Pressure>& P,
		      const boost::shared_ptr<Temperature>& T,
		      const DoubleWithUnit& Latitude, 
		      const DoubleWithUnit& Surface_height);
  virtual void add_observer(Observer<Altitude>& Obs);
  virtual void remove_observer(Observer<Altitude>& Obs);
  virtual AutoDerivativeWithUnit<double> 
  altitude(const AutoDerivativeWithUnit<double>& P) const;
  virtual AutoDerivativeWithUnit<double> 
  gravity(const AutoDerivativeWithUnit<double>& P) const;
  boost::shared_ptr<Altitude> clone() const;
  %pickle_serialization();
};

}




