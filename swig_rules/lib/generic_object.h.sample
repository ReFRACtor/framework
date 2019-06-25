#ifndef GEOCAL_GENERIC_OBJECT_H
#define GEOCAL_GENERIC_OBJECT_H
// These get defined in geocal_config.h. Undef them to stop warnings
// if we change these.
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "geocal_config.h"

namespace boost {
  namespace serialization {
    class access;
  }
}

#ifdef GEOCAL_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/export.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/version.hpp>
#define GEOCAL_EXPORT_KEY(NAME) BOOST_CLASS_EXPORT_KEY(GeoCal::NAME)
#define GEOCAL_EXPORT_KEY2(NAME, KEY) BOOST_CLASS_EXPORT_KEY2(GeoCal::NAME, KEY)
#define GEOCAL_SPLIT_MEMBER() BOOST_SERIALIZATION_SPLIT_MEMBER()
#define GEOCAL_DONT_TRACK(NAME) BOOST_CLASS_TRACKING(GeoCal::NAME, boost::serialization::track_never)
#define GEOCAL_CLASS_VERSION(NAME, N) BOOST_CLASS_VERSION(GeoCal::NAME, N)
#else
#define GEOCAL_EXPORT_KEY(NAME) /* Noop */
#define GEOCAL_EXPORT_KEY2(NAME,KEY) /* Noop */
#define GEOCAL_SPLIT_MEMBER() /* Noop */
#define GEOCAL_DONT_TRACK(NAME) /* Noop */
#define GEOCAL_CLASS_VERSION(NAME, N) /* Noop */
#endif


namespace GeoCal {
/****************************************************************//**
  For use with SWIG, it is useful to have a base class that 
  everything can be cast to. This class doesn't provide any
  functionality, other than allowing casts.
*******************************************************************/
class GenericObject {
public:
  // Have a virtual member function, which forces RTTI information to
  // be available.
  virtual ~GenericObject() {}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

GEOCAL_EXPORT_KEY(GenericObject);

#endif
