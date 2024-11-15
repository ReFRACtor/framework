#include "aerosol_shape_gaussian.h"
#include "fp_serialize_support.h"
#include "state_mapping_gaussian.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void AerosolShapeGaussian::serialize(Archive & ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AerosolExtinctionLevel);
}

FP_IMPLEMENT(AerosolShapeGaussian);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AerosolShapeGaussian, AerosolExtinction)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
     const blitz::Array<double, 1>&,
     const std::string&,
     const bool&>())
REGISTER_LUA_END()
#endif

// See base class for description

AerosolShapeGaussian::AerosolShapeGaussian(const boost::shared_ptr<Pressure>& Press,
             const blitz::Array<double, 1>& Coeffs,
             const std::string& Aerosol_name,
             const bool Linear_AOD)
  : AerosolExtinctionLevel(Press, Coeffs, Aerosol_name, boost::make_shared<StateMappingGaussian>(Press, Linear_AOD)) {}

// See base class for description
boost::shared_ptr<AerosolExtinction> AerosolShapeGaussian::clone() const
{
    boost::shared_ptr<StateMappingGaussian> gaussian_map = boost::static_pointer_cast<StateMappingGaussian>(mapping);
    return boost::shared_ptr<AerosolExtinction>
    (new AerosolShapeGaussian(press->clone(), coeff.value(),
                              aerosol_name(), gaussian_map->is_linear_total()));
}

