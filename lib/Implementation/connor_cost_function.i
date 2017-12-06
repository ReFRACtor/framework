// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "connor_cost_function.h"
%}
%base_import(cost_function)
%import "forward_model.i"
%import "state_vector.i"
%import "instrument_measurement.i"
%fp_shared_ptr(FullPhysics::ConnorCostFunction);

namespace FullPhysics {
class ConnorCostFunction : public CostFunction {
public:
  ConnorCostFunction(
          const boost::shared_ptr<StateVector>& Sv,
          const boost::shared_ptr<ForwardModel>& fm,
          const boost::shared_ptr<InstrumentMeasurement>& inst_meas);
  virtual void cost_function(const blitz::Array<double, 1>& X,
			     blitz::Array<double, 1>& OUTPUT,
			     blitz::Array<double, 1>& OUTPUT,
			     blitz::Array<double, 2>& OUTPUT) const;
};
}
