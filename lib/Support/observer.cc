#include "observer.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void CacheInvalidatedObserver::serialize(Archive & ar,
			const unsigned int version)
{
  FP_GENERIC_BASE(CacheInvalidatedObserver);
  std::string p = "empty";
  ar & FP_NVP2("placeholder", p);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void CacheInvalidatedObserver::save(Archive & , const unsigned int ) const
{
  add_ptr_serialized_reference(this);
}

template<class Archive> 
void CacheInvalidatedObserver::load(Archive & , const unsigned int) 
{ 
} 

FP_IMPLEMENT(CacheInvalidatedObserver);
#endif
