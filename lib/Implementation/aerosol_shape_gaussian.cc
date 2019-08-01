#include "aerosol_shape_gaussian.h"
#include "mapping_log.h"
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

const double AerosolShapeGaussian::min_aod = 1e-9;
// See base class for description

// TODO: Add tests for clone methods to verify things can be cloned without changes
//boost::shared_ptr<AerosolExtinction> AerosolShapeGaussian::clone
//(const boost::shared_ptr<Pressure>& Pres) const
//{
//    return boost::shared_ptr<AerosolExtinction>
//           (new AerosolShapeGaussian(Pres, used_flag, coeff.value(),
//                                     aerosol_name(), linear_aod));
//}

AerosolShapeGaussian::AerosolShapeGaussian(const boost::shared_ptr<Pressure>& Press,
             const blitz::Array<bool, 1>& Flag,
             const blitz::Array<double, 1>& Coeffs,
             const std::string& Aerosol_name,
             const bool Linear_AOD)
  : AerosolExtinctionLevel(Press, Flag, Coeffs, Aerosol_name, boost::make_shared<MappingGaussian>(Press)),
                           linear_aod(Linear_AOD) {}
