#ifndef BOOST_RATIONAL_SERIALIZE_SUPPORT_H
#define BOOST_RATIOANL_SERIALIZE_SUPPORT_H
#include "fp_serialize_support.h"

// boost::rational doesn't already have a serialization defined, so
// we'll go ahead and define one.
#include <boost/rational.hpp>
namespace boost {
  namespace serialization {
    template <typename Archive, typename T>
    void save(Archive& ar, const boost::rational<T>& R,
	      const unsigned int version);
    template <typename Archive, typename T>
    void load(Archive& ar, boost::rational<T>& R,
	      const unsigned int version);
  }
}

typedef boost::rational<int> boost_rational_int;

BOOST_SERIALIZATION_SPLIT_FREE(boost_rational_int);

#endif
