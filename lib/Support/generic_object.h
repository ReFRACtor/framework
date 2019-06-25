#ifndef GENERIC_OBJECT_H
#define GENERIC_OBJECT_H
#include "fp_config.h"

namespace boost {
  namespace serialization {
    class access;
  }
}

#ifdef FP_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/export.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/version.hpp>
#define FP_EXPORT_KEY(NAME) BOOST_CLASS_EXPORT_KEY(FullPhysics::NAME)
#define FP_EXPORT_KEY2(NAME, KEY) BOOST_CLASS_EXPORT_KEY2(FullPhysics::NAME, KEY)
#define FP_SPLIT_MEMBER() BOOST_SERIALIZATION_SPLIT_MEMBER()
#define FP_DONT_TRACK(NAME) BOOST_CLASS_TRACKING(FullPhysics::NAME, boost::serialization::track_never)
#define FP_CLASS_VERSION(NAME, N) BOOST_CLASS_VERSION(FullPhysics::NAME, N)
#else
#define FP_EXPORT_KEY(NAME) /* Noop */
#define FP_EXPORT_KEY2(NAME,KEY) /* Noop */
#define FP_SPLIT_MEMBER() /* Noop */
#define FP_DONT_TRACK(NAME) /* Noop */
#define FP_CLASS_VERSION(NAME, N) /* Noop */
#endif

// Mark a function or variable as unused. This isn't usually all that
// important, but it can be useful to silence these warnings so we get
// messages about variables that are unused that we don't expect.
#ifdef __GNUC__
#  define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#  define UNUSED(x) UNUSED_ ## x
#endif

#ifdef __GNUC__
#  define UNUSED_FUNCTION(x) __attribute__((__unused__)) UNUSED_ ## x
#else
#  define UNUSED_FUNCTION(x) UNUSED_ ## x
#endif

namespace FullPhysics {
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

FP_EXPORT_KEY(GenericObject);
#endif
