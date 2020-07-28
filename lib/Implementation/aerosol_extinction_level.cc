#include "aerosol_extinction_level.h"
#include "fp_serialize_support.h"
#include "fp_exception.h"
#include "ostream_pad.h"
#include <boost/lexical_cast.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AerosolExtinctionLevel::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AerosolExtinctionImpBase);
}

FP_IMPLEMENT(AerosolExtinctionLevel);
#endif

// See base class for description
boost::shared_ptr<AerosolExtinction> AerosolExtinctionLevel::clone() const
{
  return boost::shared_ptr<AerosolExtinction>
    (new AerosolExtinctionLevel(press->clone(), used_flag, coeff.value(),
                                aerosol_name(), mapping->clone()));
}

void AerosolExtinctionLevel::calc_aerosol_extinction() const
{
  aext = mapping->fm_view(coeff);
}

void AerosolExtinctionLevel::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "    ");
  Os << "AerosolExtinctionLevel" + mapping->name() + ":\n"
     << "  Aerosol name:       " << aerosol_name() << "\n"
     << "  Coefficient:\n";
  opad << coeff.value() << "\n";
  opad.strict_sync();
  Os << "  Retrieval flag:\n";
  opad << used_flag << "\n";
  opad.strict_sync();
}
