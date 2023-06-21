#ifndef OBSERVATION_H
#define OBSERVATION_H

#include "printable.h"
#include "stacked_radiance_mixin.h"

namespace FullPhysics {
/****************************************************************//**
  This is the observed data, e.g. the L1B measured radiance data from
  the OCO-2 instrument. This can be compared against the radiance
  from a ForwardModel during a retrieval.
 *******************************************************************/

class Observation : public StackedRadianceMixin {
public:
  /// Number of spectral channels
  virtual int num_channels() const = 0;

  /// The spectral grid of the radiance values, implemented by inheriting class
  virtual SpectralDomain spectral_domain(int sensor_index) const = 0;

  /// Measured spectrum for the given spectral channel. Note that this may be empty.
  ///
  /// \param sensor_index Band to give value for
  /// \param skip_jacobian If true, don't do the Jacobian
  /// \return The set of radiances, possibly empty.
  virtual Spectrum radiance(int sensor_index, bool skip_jacobian = false) const = 0;

  virtual void print(std::ostream& Os) const
  {
    Os << "Observation";
  }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(Observation);
#endif
