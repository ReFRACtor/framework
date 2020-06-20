#ifndef OBSERVER_SERIALIZE_SUPPORT_H
#define OBSERVER_SERIALIZE_SUPPORT_H
#include "observer.h"
#include "fp_serialize_support.h"
#include <set>

namespace FullPhysics {
/****************************************************************//**
  Small class that just saves a copy of all the pointers as we 
  do serialization. We have this as a singleton for now, I'm pretty
  sure we would never want multiple versions.
*******************************************************************/

class ObserverSerializedMarker {
public:
  static ObserverSerializedMarker& instance();
  template<class T> static void add_serialized_reference(const Observer<T>* P)
  {
    instance().data.insert(static_cast<const GenericObject*>(P));
  }
  template<class T> static bool is_serialized(const Observer<T>* P)
  {
    return instance()._is_serialized(P);
  }
  static void clear()
  {
    instance().data.clear();
  }
private:
  template<class T> bool _is_serialized(const Observer<T>* P)
  { return data.find(static_cast<const GenericObject*>(P)) != data.end(); }
  static ObserverSerializedMarker* instance_;
  ObserverSerializedMarker() {}
  std::set<const GenericObject*> data;
};
}

#define FP_OBSERVER_SERIALIZE(NAME) \
template<> template<class Archive>	      \
void Observer ## NAME::serialize(Archive& ar, \
			 const unsigned int version) \
{ \
  FP_GENERIC_BASE(Observer ## NAME); \
  std::string p = "empty"; \
  ar & FP_NVP2("placeholder", p); \
  boost::serialization::split_member(ar, *this, version); \
} \
template<> template<class Archive> \
void Observer ## NAME::save(Archive & , const unsigned int ) const \
{ \
  ObserverSerializedMarker::add_serialized_reference(this); \
} \
template<> template<class Archive> \
void Observer ## NAME::load(Archive & ar, const unsigned int)	\
{ \
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
void Observable ## NAME::save(Archive & ar, const unsigned int ) const \
{ \
  std::list<boost::weak_ptr<Observer ## NAME > > olist_keep; \
  BOOST_FOREACH(const boost::weak_ptr<Observer ## NAME >& t, olist) { \
    boost::shared_ptr<Observer ## NAME> t2 = t.lock(); \
    if(t2 && ObserverSerializedMarker::is_serialized(t2.get())) \
      olist_keep.push_back(t); \
  } \
  std::cerr << "olist_keep.size() " << olist_keep.size() << "\n"; \
  ar & FP_NVP2("olist", olist_keep); \
  std::cerr << "done with olist\n"; \
} \
template<> template<class Archive> \
void Observable ## NAME::load(Archive & ar, const unsigned int)	\
{ \
  ar & FP_NVP(olist); \
} \
FP_IMPLEMENT(Observer ## NAME); \
FP_IMPLEMENT(Observable ## NAME);

#endif
