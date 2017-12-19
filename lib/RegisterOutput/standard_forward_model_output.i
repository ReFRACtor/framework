// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "standard_forward_model_output.h"
#include "ils_instrument.h"
#include "pressure.h"

// Needed for type conversions in SWIG
#include "sub_state_vector_array.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "standard_forward_model.i"

%fp_shared_ptr(FullPhysics::StandardForwardModelOutput);

namespace FullPhysics {
class StandardForwardModelOutput : public RegisterOutputBase {
public:
  StandardForwardModelOutput(const boost::shared_ptr<StandardForwardModel>& Fm, const boost::shared_ptr<Observation>& inst_meas);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
};
}
