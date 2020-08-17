#include "aerosol_extinction_linear.h"
#include "fp_serialize_support.h"
#include "state_mapping_linear.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AerosolExtinctionLinear::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AerosolExtinctionLevel);
}

FP_IMPLEMENT(AerosolExtinctionLinear);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AerosolExtinctionLinear, AerosolExtinction)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
                          const blitz::Array<double, 1>&,
                          const std::string&>())
REGISTER_LUA_END()
#endif

AerosolExtinctionLinear::AerosolExtinctionLinear(const boost::shared_ptr<Pressure>& Press,
              const blitz::Array<double, 1>& Aext,
              const std::string& Aerosol_name)
    : AerosolExtinctionLevel(Press, Aext, Aerosol_name, boost::make_shared<StateMappingLinear>()) {}

// See base class for description
boost::shared_ptr<AerosolExtinction> AerosolExtinctionLinear::clone() const
{
    return boost::shared_ptr<AerosolExtinction>
    (new AerosolExtinctionLinear(press->clone(), coeff.value(), aerosol_name()));
}

