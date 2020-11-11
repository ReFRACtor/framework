#include "temperature_met.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void TemperatureMet::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TemperatureOffset)
    & FP_NVP(met);
}

FP_IMPLEMENT(TemperatureMet);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(TemperatureMet, Temperature)
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
			  const boost::shared_ptr<Pressure>&,
			  double>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Create an Temperature. 
//-----------------------------------------------------------------------

TemperatureMet::TemperatureMet
(const boost::shared_ptr<Meteorology>& Met_file,
 const boost::shared_ptr<Pressure>& Press,
 double Temp_offset)
: TemperatureOffset(Press, Temp_offset), met(Met_file)
{
}

// See base class for description of this function
boost::shared_ptr<Temperature> 
TemperatureMet::clone() const
{
  boost::shared_ptr<Temperature> res
    (new TemperatureMet(met, press->clone(), coefficient()(0).value()));
  return res;
}

void TemperatureMet::print(std::ostream& Os) const 
{
  OstreamPad opad(Os, "    ");
  Os << "TemperatureMet:\n"
     << "  Temperature offset: " << temperature_offset() << "\n"
     << "  Meteorology:\n";
  opad << *met << "\n";
  opad.strict_sync();
}
