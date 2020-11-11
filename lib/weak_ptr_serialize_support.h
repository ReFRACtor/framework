#ifndef WEAK_PTR_SERIALIZED_SUPPORT_H
#define WEAK_PTR_SERIALIZED_SUPPORT_H
#include "generic_object.h"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_free.hpp>

/****************************************************************//**
  There is a complication in serializing a boost::weak_ptr. Normally
  when you encounter a boost::shared_ptr you want to serialize the
  object, e.g. object A contains a boost:shared_ptr pointing to object
  B. 

  boost::weak_ptr on the other hand is a looser connection. It means
  an object A has a reference to access B, but doesn't really "own"
  it (e.g., B is a Observer of A, and A notifies it about changes).

  There is no one "right" way to handle this, it really just depends
  on what the semantics are for you objects.

  There are two ways to handled this already in boost serialization:

  1. You don't want to keep any references. The weak_ptr is just not
     included in the serialization of class A.
  2. You always want the weak_ptr included just like we do with
     boost::shared_ptr. boost has serialization code that does this
     in boost/serialization/weak_ptr.hpp.

  Often neither of these ways are what you want. For example, in 
  ReFRACtor we have Pressure, which is Observerable. For case (1)
  we loose all internal connections in something like a ForwardModel,
  objects which need to know then Pressure changes loose that
  connection when we serialize the ForwardModel. On the other hand,
  we don't want to serialize just a Pressure object and also get the 
  entire world through a series of pointers. Also there might be
  observers like a logger writing out when Pressure changes that we
  wouldn't want to serialize.

  We supply an alternative to the two where we include all weak_ptr
  to objects that are otherwise included in a serialization. So 
  from our example, serializing a single Pressure object would have 
  no weak_ptr serialized, while a ForwardModel would have all weak_ptr
  from Pressure to other objects in the ForwardModel, but not to
  object not in ForwardModel.

  This is accomplished by doing the serialization twice. The first 
  time, we mark all the GenericObject  that get serialized. The second
  time, we keep all weak_ptr that point to one of these objects.

  The serialization code needs to "know" to call the 
  add_serialized_reference function in this header. This isn't too
  hard in practice, usually only special objects can have weak_ptr to
  them, e.g Observers of a Observable. There is sample code for this
  in ReFRACtor if you want to look at an example of its use.

  Note that the GenericObjectMap can be used for collections of
  objects you want to have serialized together that don't already have
  a containing object, e.g. a ForwardModel plus a StateVector.

  Also note that you are not required to use this particular method if it
  doesn't fit your semantics - you can use one of the others supplied
  by boost or create your own way to handle this.
*******************************************************************/

namespace SWIG_MAPPER_NAMESPACE {
  void clear_ptr_serialized_reference();
  void add_ptr_serialized_reference(const GenericObject* P);
  bool is_ptr_serialized(const GenericObject* P);  
}

namespace boost {
namespace serialization{

template<class Archive, class T>
inline void save
(Archive & ar,
 const boost::weak_ptr< T > &t,
 const unsigned int /* file_version */
)
{
  boost::shared_ptr< T > sp = t.lock();
  const SWIG_MAPPER_NAMESPACE::GenericObject* ptr = 0;
  if(sp) {
    ptr = dynamic_cast<const SWIG_MAPPER_NAMESPACE::GenericObject*>(sp.get());
    if(ptr && !SWIG_MAPPER_NAMESPACE::is_ptr_serialized(ptr))
      sp.reset();
  }
  ar << boost::serialization::make_nvp("weak_ptr", sp);
}

template<class Archive, class T>
inline void load
(Archive & ar,
 boost::weak_ptr< T > &t,
 const unsigned int /* file_version */
)
{
  boost::shared_ptr< T > sp;
  ar >> boost::serialization::make_nvp("weak_ptr", sp);
  t = sp;
}

template<class Archive, class T>
inline void serialize
(
 Archive & ar,
 boost::weak_ptr< T > &t,
 const unsigned int file_version
 )
{
  boost::serialization::split_free(ar, t, file_version);
}
}
}

#endif
