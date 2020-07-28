#include "atmosphere_fixture.h"
using namespace FullPhysics;
using namespace blitz;

AtmosphereFixture::AtmosphereFixture()
{
  atm = dynamic_cast<const AtmosphereStandard&>(*config_atmosphere).clone();

  // Create a new SV for our cloned Atmosphere to use
  statev.reset(new StateVector());

  attach_atmosphere_to_sv(atm, statev);
}

void AtmosphereFixture::attach_atmosphere_to_sv(boost::shared_ptr<AtmosphereStandard>& atmosphere, boost::shared_ptr<StateVector>& state_vector)
{
  // Attach state vector to atmosphere
  state_vector->add_observer(*atmosphere);

  // Need to reattach the cloned atmosphere components to the SV one by one as done in the Lua config and in the same order
  atmosphere->attach_children_to_sv(*state_vector);

  // Push values initial values into the state vector
  state_vector->update_state(config_initial_guess->initial_guess());
}
