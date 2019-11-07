#include <boost/make_shared.hpp>
#include "absorber_vmr_level_scaled.h"
#include "mapping_scale.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrLevelScaled, AbsorberVmr)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
			  const blitz::Array<double, 1>&,
			  double, 
			  bool,
			  const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
AbsorberVmrLevelScaled::AbsorberVmrLevelScaled(const boost::shared_ptr<Pressure>& Press,
 const blitz::Array<double, 1>& Vmr_profile,
 double Scale,                         
 bool Scale_flag,
 const std::string& Gas_name)
: AbsorberVmrLevel(Press, Vmr_profile, Scale_flag, Gas_name, boost::make_shared<MappingScale>(Scale, Vmr_profile.copy()))
{ 
}

double AbsorberVmrLevelScaled::scale_factor() const { return coeff(0).value(); }

boost::shared_ptr<AbsorberVmr> AbsorberVmrLevelScaled::clone() const
{
    return boost::shared_ptr<AbsorberVmr>
    (new AbsorberVmrLevelScaled(press->clone(), vmr_profile(), coeff(0).value(),used_flag(0),
                                gas_name()));
}
