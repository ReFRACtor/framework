#include "aerosol_extinction_log.h"
#include "fp_serialize_support.h"
#include "state_mapping_log.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AerosolExtinctionLog::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AerosolExtinctionLevel);
}

FP_IMPLEMENT(AerosolExtinctionLog);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AerosolExtinctionLog, AerosolExtinction)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
                          const blitz::Array<bool, 1>&, 
                          const blitz::Array<double, 1>&,
                          const std::string&>())
REGISTER_LUA_END()
#endif


AerosolExtinctionLog::AerosolExtinctionLog(const boost::shared_ptr<Pressure>& Press,
                                           const blitz::Array<bool, 1>& Flag,
                                           const blitz::Array<double, 1>& Aext,
                                           const std::string& Aerosol_name)
    : AerosolExtinctionLevel(Press, Flag, Aext, Aerosol_name, boost::make_shared<StateMappingLog>()) {}

// See base class for description
boost::shared_ptr<AerosolExtinction> AerosolExtinctionLog::clone() const
{
    return boost::shared_ptr<AerosolExtinction>
    (new AerosolExtinctionLog(press->clone(), used_flag, blitz::Array<double, 1>(exp(coeff.value())), aerosol_name()));
}

