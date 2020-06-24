#include "pressure_level_input.h"
#include "ostream_pad.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void PressureLevelInput::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  FP_GENERIC_BASE(PressureLevelInput);
  ar & FP_NVP(press_level);
}

FP_IMPLEMENT(PressureLevelInput);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(PressureLevelInput)
.def(luabind::constructor<const blitz::Array<double, 1>&>())
.def(luabind::constructor<const HdfFile&>())
.def(luabind::constructor<const HdfFile&, const std::string>())
.def("pressure_level", &PressureLevelInput::pressure_level)
REGISTER_LUA_END()
#endif

void PressureLevelInput::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "  ");
  Os << "PressureLevelInput:\n";
  opad << press_level;
  opad.strict_sync();
}
