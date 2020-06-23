#include "generic_object_map.h"
#include "geocal_serialize_support.h"

using namespace GeoCal;

template<class Archive>
void GenericObjectMap::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  GEOCAL_GENERIC_BASE(GenericObjectMap);
  typedef std::map<std::string, boost::shared_ptr<GenericObject> >
    MapGenericObject;
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MapGenericObject);
}

GEOCAL_IMPLEMENT(GenericObjectMap);
