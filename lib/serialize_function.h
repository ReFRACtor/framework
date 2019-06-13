#ifndef SERIALIZE_FUNCTION_H
#define SERIALIZE_FUNCTION_H
#include "generic_object.h"
#include <string>
#include <boost/shared_ptr.hpp>

namespace SWIG_MAPPER_NAMESPACE {
bool have_serialize_supported();
void serialize_write(const std::string& Fname, 
		     const boost::shared_ptr<GenericObject>& Obj);
void serialize_write_binary(const std::string& Fname, 
		     const boost::shared_ptr<GenericObject>& Obj);

std::string serialize_write_string(const boost::shared_ptr<GenericObject>& Obj);


boost::shared_ptr<GenericObject> 
serialize_read_generic(const std::string& Fname);
boost::shared_ptr<GenericObject> 
serialize_read_binary_generic(const std::string& Fname);

boost::shared_ptr<GenericObject> 
serialize_read_generic_string(const std::string& Data);

std::string serialize_write_binary(const boost::shared_ptr<GenericObject>& Obj);
boost::shared_ptr<GenericObject> 
serialize_read_binary(const std::string& Data);

template<class T> inline boost::shared_ptr<T> 
serialize_read(const std::string& Fname)
{
  return boost::dynamic_pointer_cast<T>(serialize_read_generic(Fname));
}

template<class T> inline boost::shared_ptr<T> 
serialize_read_binary(const std::string& Fname)
{
  return boost::dynamic_pointer_cast<T>(serialize_read_binary_generic(Fname));
}
  
template<class T> inline boost::shared_ptr<T> 
serialize_read_string(const std::string& Data)
{
  return boost::dynamic_pointer_cast<T>(serialize_read_generic_string(Data));
}
}
#endif
