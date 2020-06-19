#ifndef OBSERVER_SERIALIZE_SUPPORT_H
#define OBSERVER_SERIALIZE_SUPPORT_H
#include "observer.h"
#include "fp_serialize_support.h"

#define FP_OBSERVER_SERIALIZE(NAME) \
template<> template<class Archive>	      \
void Observer ## NAME::serialize(Archive& ar, \
			 const unsigned int UNUSED(version)) \
{ \
  FP_GENERIC_BASE(Observer ## NAME); \
  std::string p = "empty"; \
  ar & FP_NVP2("placeholder", p); \
} \
template<> template<class Archive> \
void Observable ## NAME::serialize(Archive& ar, \
			      const unsigned int version) \
{ \
  FP_GENERIC_BASE(Observable ##NAME); \
  ar & FP_NVP(ref_list); \
  boost::serialization::split_member(ar, *this, version); \
} \
template<> template<class Archive> \
void Observable ## NAME::save(Archive & , const unsigned int ) const \
{ \
} \
template<> template<class Archive> \
void Observable ## NAME::load(Archive & , const unsigned int)	\
{ \
    BOOST_FOREACH(boost::shared_ptr<Observer ## NAME >& t, ref_list) { \
      add_observer(*t); \
    } \
} \
FP_IMPLEMENT(Observer ## NAME); \
FP_IMPLEMENT(Observable ## NAME);

#endif
