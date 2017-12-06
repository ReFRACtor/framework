// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "spectral_parameters_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "forward_model.i"

%fp_shared_ptr(FullPhysics::SpectralParametersOutput);

namespace FullPhysics {
class SpectralParametersOutput : public RegisterOutputBase {
public:
  SpectralParametersOutput(const boost::shared_ptr<ForwardModel>& Fm);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
};
}


