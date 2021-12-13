#ifndef ATMOSPHERE_FIXTURE_H
#define ATMOSPHERE_FIXTURE_H
#include "lua_configuration_fixture.h"
#include "atmosphere_standard.h"
#include "state_vector.h"

namespace FullPhysics {
/****************************************************************//**
  This is a test fixture that creates an Atmosphere and related
  Statevector based on the unit test data set up in 
  ConfigurationFixture.

  This is pretty similar to ConfigurationFixture, but we have created 
  a state vector that just has a clone of the Atmosphere. The
  Atmsosphere clone gets created each time this Fixture is created, so
  we have a clean Atmosphere for testing (a previous test might have
  purposely put the Atmosphere in an error state). If you don't need a
  clean new Atmosphere each time, you might just want to use
  ConfigurationFixture.
*******************************************************************/
class AtmosphereFixture : public ConfigurationFixture {
public:
  AtmosphereFixture();
  virtual ~AtmosphereFixture() {statev->remove_observer(*atm);}

  /// Atmosphere read from ConfigurationFixture
  boost::shared_ptr<AtmosphereStandard> atm;

  /// Helper routine for initializing an Atmosphere in such a way that jacobians are set up and linked to a StateVector
  void attach_atmosphere_to_sv(boost::shared_ptr<AtmosphereStandard>& atmosphere, boost::shared_ptr<StateVector>& state_vector);

  void set_surface_pressure(double x);
 
  /// Statevector that atm is attached as an observer
  boost::shared_ptr<StateVector> statev;

  int state_vector_size() const 
  {return statev->state().extent(blitz::firstDim);}
};

}
#endif
