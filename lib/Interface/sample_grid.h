#ifndef SAMPLE_GRID_H
#define SAMPLE_GRID_H

#include "state_vector_observer.h"
#include "observer.h"
#include "spectral_domain.h"

namespace FullPhysics {
/****************************************************************//**
  This class calculates the wavenumber for each sample in a single
  band of an Instrument.

*******************************************************************/
class SampleGrid : virtual public StateVectorObserver,
              public Observable<SampleGrid> {

public:
  virtual ~SampleGrid() {}
  virtual void add_observer(Observer<SampleGrid>& Obs)
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<SampleGrid>& Obs)
  { remove_observer_do(Obs, *this);}

  //-----------------------------------------------------------------------
  /// Clone a SampleGrid object. Note that the cloned version will *not*
  /// be attached to and StateVector or Observer<SampleGrid>, although you
  /// can of course attach them after receiving the cloned object.
  ///
  /// Because this isn't attached to the StateVector, one use of the
  /// clone operator is to create a "frozen" SampleGrid object.
  //-----------------------------------------------------------------------

    virtual boost::shared_ptr<SampleGrid> clone() const = 0;

  //-----------------------------------------------------------------------
  /// Returns as list of grid points for each instrument sample, and the
  /// gradient of the points wrt the state vector. This is for the
  /// full instrument samples, i.e., any windowing etc. happens in later
  /// processing.
  //-----------------------------------------------------------------------

  virtual SpectralDomain sample_grid() const = 0;

  //-----------------------------------------------------------------------
  /// Supports old Dispersion call interface
  /// gradient of the points wrt the state vector. This is for the
  /// full instrument samples, i.e., any windowing etc. happens in later
  /// processing.
  //-----------------------------------------------------------------------
  virtual SpectralDomain pixel_grid() const { return this->sample_grid(); }
};
}
#endif
