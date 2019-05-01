#include <boost/make_shared.hpp>
#include "absorber_vmr_level_scaled2.h"
#include "mapping_scale.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------
AbsorberVmrLevelScaled2::AbsorberVmrLevelScaled2(const boost::shared_ptr<Pressure>& Press,
                                                 const blitz::Array<double, 1>& Vmr_profile,
                                                 double Scale,
                                                 const blitz::Array<bool, 1>& Scale_flag,
                                                 const std::string& Gas_name)
: AbsorberVmrLevel(Press, Vmr_profile, Scale_flag, Gas_name, boost::make_shared<MappingScale>(Scale, Vmr_profile.copy()))
{
}

