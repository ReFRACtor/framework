#include <boost/make_shared.hpp>
#include "absorber_vmr_level_scaled.h"
#include "fp_serialize_support.h"
#include "state_mapping_scale.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AbsorberVmrLevelScaled::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbsorberVmrLevel);
}

FP_IMPLEMENT(AbsorberVmrLevelScaled);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrLevelScaled, AbsorberVmr)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
                          const blitz::Array<double, 1>&,
                          double, 
                          const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
AbsorberVmrLevelScaled::AbsorberVmrLevelScaled(const boost::shared_ptr<Pressure>& Press,
 const blitz::Array<double, 1>& Vmr_profile,
 double Scale,                         
 const std::string& Gas_name)
: AbsorberVmrLevel(Press, Vmr_profile, Gas_name, boost::make_shared<StateMappingScale>(Scale, Vmr_profile.copy()))
{ 
}

double AbsorberVmrLevelScaled::scale_factor() const { return coeff(0).value(); }

boost::shared_ptr<AbsorberVmr> AbsorberVmrLevelScaled::clone() const
{
    return boost::shared_ptr<AbsorberVmr>
    (new AbsorberVmrLevelScaled(mapped_pressure->clone(), vmr_profile(), coeff(0).value(), gas_name()));
}

std::string AbsorberVmrLevelScaled::state_vector_name_i(int coeff_idx) const
{
  range_check(coeff_idx, 0, 1);
  std::stringstream sv_name;
  sv_name << gas_name() << " Scaled VMR";
  return sv_name.str();
}
