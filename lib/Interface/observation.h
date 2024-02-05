#ifndef OBSERVATION_H
#define OBSERVATION_H

#include "printable.h"
#include "stacked_radiance_mixin.h"
#include "observer.h"

namespace FullPhysics {
/****************************************************************//**
  This is the observed data, e.g. the L1B measured radiance data from
  the OCO-2 instrument. This can be compared against the radiance
  from a ForwardModel during a retrieval.

  Note that the interface is that this class returns the set of
  observations that we are fitting for - whatever is needed in
  ModelMeasure. So for example, this might by only for a subset of the
  values in file, and might exclude bad pixels (e.g., like using a
  SpectralWindow to select a subset).

  We had considered doing the subsetting of the data outside of this
  class, but there seemed to be too many variations to this (e.g., is
  the underlying file already subsetted for spectral range?  Does it
  have bad pixels in it, or are these filtered out?). So it is the
  responsibility of this class to handle that (e.g, it might take a
  SpectralWindow). For the StandardForwardModel, this should match the
  low_resolution_grid from the ForwardModelSpectralGrid, and in
  general it should match the spectral domain of the model in a
  ModelMeasure.
 *******************************************************************/

class Observation : public StackedRadianceMixin,
		    public Observable<Observation> {
public:
  virtual ~Observation() {}
  
  virtual void add_observer(Observer<Observation>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Observation>& Obs) 
  { remove_observer_do(Obs, *this);}

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
FP_CLASS_VERSION(Observation, 1);
FP_EXPORT_OBSERVER_KEY(Observation);
#endif
