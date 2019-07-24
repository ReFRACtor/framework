#include "aerosol_extinction_log.h"
#include "mapping_log.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

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
    : AerosolExtinctionLevel(Press, Flag, Aext, Aerosol_name, boost::make_shared<MappingLog>()) {}
