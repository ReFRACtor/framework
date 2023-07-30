// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "observation_sv.h"
%}

%base_import(state_vector)
%base_import(state_vector_observer)
%import "spectrum.i"

%fp_shared_ptr(FullPhysics::ObservationSv)

namespace FullPhysics {
class ObservationSv : public Observation,
		      public StateVectorObserver {
public:
  virtual SpectralDomain spectral_domain(int sensor_index) const = 0;
  virtual boost::shared_ptr<ObservationSv> clone() const;
  %pickle_serialization();
};
}

