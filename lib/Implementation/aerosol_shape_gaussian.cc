#include "aerosol_shape_gaussian.h"
#include "mapping_gaussian.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AerosolShapeGaussian, AerosolExtinction)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
     const blitz::Array<bool, 1>&,
     const blitz::Array<double, 1>&,
     const std::string&,
     const bool&>())
REGISTER_LUA_END()
#endif

// See base class for description

// TODO: Add tests for clone methods to verify things can be cloned without changes
boost::shared_ptr<AerosolExtinction> AerosolShapeGaussian::clone
(const boost::shared_ptr<Pressure>& Pres) const
{
//    This works.
//    boost::shared_ptr<MappingGaussian> mapping_gaussian = boost::static_pointer_cast<MappingGaussian>(mapping);
//    return boost::shared_ptr<AerosolExtinction>
//           (new AerosolShapeGaussian(Pres, used_flag, coeff.value(),
//                                     aerosol_name(), mapping_gaussian->is_linear_total() ));
//

//  TODO: This does not work -- default implementation from AerosolExtinctionLevel
    return boost::shared_ptr<AerosolExtinction>
      (new AerosolExtinctionLevel(Pres, used_flag, coeff.value(), aerosol_name(), mapping));
}

AerosolShapeGaussian::AerosolShapeGaussian(const boost::shared_ptr<Pressure>& Press,
             const blitz::Array<bool, 1>& Flag,
             const blitz::Array<double, 1>& Coeffs,
             const std::string& Aerosol_name,
             const bool Linear_AOD)
  : AerosolExtinctionLevel(Press, Flag, Coeffs, Aerosol_name, boost::make_shared<MappingGaussian>(Press, Linear_AOD)) {}
