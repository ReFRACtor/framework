#ifndef GENERIC_OBJECT_MAP_H
#define GENERIC_OBJECT_MAP_H
#include "printable.h"
#include <boost/shared_ptr.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/assign.hpp>
#include <algorithm>
#include <string>
#include <map>

namespace GeoCal {

/****************************************************************//**
  Boost serialization only maintains pointers for a single object
  being serialized. So for example we may have object 'a' and object
  'b'. Object 'b' contains a pointer to object 'a'. If we serialize
  object 'b' and object 'a' separately, and then read them back in,
  we end up with two copies of 'a' - one from the serialization of the
  original 'a' and one from 'b'. If instead, we have an object c that
  contains both 'a' and 'b', then there will be only one 'a', since
  boost in serializing c realized there are two pointers that point
  to the same object.

  As a help when we don't already have a object 'c', this map contains
  any arbitrary number of objects, allowing us maintain the
  relationships.

  This serves much the same purpose as "shelve" does for "pickle" in
  python.
*******************************************************************/
class GenericObjectMap : public Printable<GenericObjectMap>,
			 public std::map<std::string,
					 boost::shared_ptr<GenericObject> >{
public:
  GenericObjectMap() {}
  virtual ~GenericObjectMap() {}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const {Os << "GenericObjectMap";}

//-----------------------------------------------------------------------
/// List of keys, for python
//-----------------------------------------------------------------------

  std::vector<std::string> keys() const
  {
    std::vector<std::string> r;
    boost::copy(*this | boost::adaptors::map_keys, std::back_inserter(r));
    return r;
  }
  
//-----------------------------------------------------------------------
/// Save object - this is something more easy to call from python.
//-----------------------------------------------------------------------
  
  void set_generic(const std::string& Key,
		   const boost::shared_ptr<GenericObject>& V)
  { (*this)[Key] = V; }
  
//-----------------------------------------------------------------------
/// Return object if we have one, or a null pointer if there is no
/// value in the map.
//-----------------------------------------------------------------------
  
  boost::shared_ptr<GenericObject> get_generic(const std::string& Key) const
  {
    std::map<std::string, boost::shared_ptr<GenericObject> >::const_iterator
      i = find(Key);
    if(i != end())
      return i->second;
    return boost::shared_ptr<GenericObject>();
  }

//-----------------------------------------------------------------------
/// Return object if we have one, or a null pointer if there is no
/// value in the map or it doesn't match type
//-----------------------------------------------------------------------
  
  template<class T> inline boost::shared_ptr<T> get(const std::string& Key)
    const
  {
    return boost::dynamic_pointer_cast<T>(get_generic(Key));
  }

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

GEOCAL_EXPORT_KEY(GenericObjectMap);
#endif
