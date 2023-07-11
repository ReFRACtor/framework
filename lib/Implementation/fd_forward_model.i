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
  %python_attribute(num_channels, int);
  virtual SpectralDomain spectral_domain(int sensor_index) const;
  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false) const;
  virtual void setup_grid();
  virtual SpectralDomain::TypePreference spectral_domain_type_preference() const;
  %python_attribute(state_vector, boost::shared_ptr<StateVector>);
  %pickle_serialization();
};
}
