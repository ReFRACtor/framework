#ifndef AEROSOL_H
#define AEROSOL_H
#include "state_vector_observer.h"
#include "accumulated_timer.h"
#include "pressure.h"
#include "relative_humidity.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the aerosol portion of the state.

  Other objects may depend on the aerosol, and should be updated
  when the aerosol is updated. To facilitate that, this class in
  an Oberverable, and objects can add themselves as Observers to be
  notified when the aerosol is updated.

  I'm not really sure what the interface for this class should be.
  Right now it is used only by AtmosphereStandard, and there is only one
  instance AerosolOptical, so the functions are what AtmosphereStandard
  needs. But we may perhaps want to modify this in the future to be
  more general. 
*******************************************************************/

class Aerosol: public StateVectorObserver, public Observable<Aerosol> {
public:
  virtual ~Aerosol() {}
  // Used as a convenience to collect timing information to report in 
  // logging
  static AccumulatedTimer timer;

  virtual void add_observer(Observer<Aerosol>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Aerosol>& Obs)
  { remove_observer_do(Obs, *this);}

  virtual boost::shared_ptr<Aerosol> clone() const = 0;

//-----------------------------------------------------------------------
/// Returns the portion of the phase function moments that come from 
/// a specific aerosol partiicle.
/// \param wn The wave number.
/// \param pidx Aerosol particle index
/// \param nummom Number of moments to fill in
/// \param numscat Number of scatters to fill in
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 3> pf_mom(double wn, int pindex, int nummom = -1, int numscat = -1) const = 0;

//-----------------------------------------------------------------------
/// This calculates the portion of the phase function moments that
/// come from the aerosol.
/// \param wn The wave number.
/// \param frac_aer This is number_active_layer() x number_particle()
/// \param nummom Number of moments to fill in
/// \param numscat Number of scatters to fill in
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 3> pf_mom(double wn, 
         const ArrayAd<double, 2>& frac_aer,
         int nummom = -1, int numscat = -1) const = 0;

//-----------------------------------------------------------------------
/// Number of aerosol particles
//-----------------------------------------------------------------------

  virtual int number_particle() const = 0;

//-----------------------------------------------------------------------
/// This gives the extinction optical depth for each layer, for the given wave
/// number. Note this only includes the aerosol portion of this,
/// Atmosphere class combines this with Absorbers and rayleigh
/// scattering.
///
/// This calculates the derivatives with respect to the state vector.
///
/// This has size of number_active_layer() x number_particle().
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 2> extinction_optical_depth_each_layer(double wn) const = 0;

//-----------------------------------------------------------------------
/// This gives the scattering optical depth for each layer, for the given wave
/// number, for the given particle. Note this only includes the
/// aerosol portion of this, 
/// Atmosphere class combines this with Rayleigh scattering.
///
/// We take in the optical depth of each layer. This is just what is
/// returned by extinction_optical_depth_each_layer(), we take this in because
/// we can change what the derivative of extinction_optical_depth_each_layer is
/// respect to, e.g. in AtmosphereStandard we use taua_i.
///
/// This calculates the derivative with respect to whatever variables
/// Od is relative to.
///
/// This has size of number_active_layer()
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 2> scattering_optical_depth_each_layer(double wn) const = 0;
};
}
#endif
