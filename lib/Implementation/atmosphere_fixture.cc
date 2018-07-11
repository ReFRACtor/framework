#include "atmosphere_fixture.h"
#include "pressure_fixed_level.h"
using namespace FullPhysics;
using namespace blitz;

AtmosphereFixture::AtmosphereFixture()
{
  atm = dynamic_cast<const AtmosphereOco&>(*config_atmosphere).clone();
  press_level = config_pressure_level_input;

  // Create a new SV for our cloned Atmosphere to use
  statev.reset(new StateVector);
  statev->add_observer(*atm);

  // Need to reattach the cloned atmosphere components to the SV one by one as done in the Lua config and in the same order
  atm->attach_children_to_sv(*statev);

  statev->update_state(config_initial_guess->initial_guess());
}

/// Set the surface pressure of the atmosphere.
void AtmosphereFixture::set_surface_pressure(double x)
{
  atm->set_surface_pressure_for_testing(x);
}

