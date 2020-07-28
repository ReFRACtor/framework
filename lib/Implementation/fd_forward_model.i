// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "fd_forward_model.h"
%}
%base_import(forward_model)
%import "state_vector.i"
%fp_shared_ptr(FullPhysics::FdForwardModel);

namespace FullPhysics {
class FdForwardModel : public ForwardModel {
public:
  FdForwardModel(const boost::shared_ptr<ForwardModel>& Real_Forward_model,
  		 const boost::shared_ptr<StateVector>& Sv,
  		 const blitz::Array<double, 1>& Perturbation);
  virtual Spectrum radiance(int channel_index, bool skip_jacobian = false) const;
  %python_attribute(state_vector, boost::shared_ptr<StateVector>);
  %pickle_serialization();
};
}
