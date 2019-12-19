#include "atmosphere_fixture.h"
#include "pressure_fixed_level.h"
using namespace FullPhysics;
using namespace blitz;

AtmosphereFixture::AtmosphereFixture()
{
  atm = dynamic_cast<const AtmosphereStandard&>(*config_atmosphere).clone();
  press_level = config_pressure_level_input;

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

/// Set the surface pressure of the atmosphere.
void AtmosphereFixture::set_surface_pressure(double x)
{
  atm->set_surface_pressure_for_testing(x);
}

