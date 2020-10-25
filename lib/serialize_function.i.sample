// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "geocal_common.i"
%{
#include "serialize_function.h"
%}
%import "generic_object.i"

namespace GeoCal {
bool have_serialize_supported();
void serialize_write(const std::string& Fname, 
		     const boost::shared_ptr<GenericObject>& Obj);
std::string serialize_write_string(const boost::shared_ptr<GenericObject>& Obj);

boost::shared_ptr<GenericObject> 
serialize_read_generic(const std::string& Fname);

boost::shared_ptr<GenericObject> 
serialize_read_generic_string(const std::string& Data);

void serialize_write_binary(const std::string& Fname, 
			    const boost::shared_ptr<GenericObject>& Obj);
boost::shared_ptr<GenericObject> 
serialize_read_binary_generic(const std::string& Fname);
}

// For binary data, we need to go to a python byte array rather than
// a python string
%typemap(out) std::string {
  $result = PyBytes_FromStringAndSize($1.data(), $1.size());
}

%typemap(in) const std::string& {
 if(PyBytes_Check($input))
     $1 = new std::string(PyBytes_AS_STRING($input), PyBytes_GET_SIZE($input));
 else
   if(PyByteArray_Check($input))
     $1 = new std::string(PyByteArray_AS_STRING($input), PyByteArray_GET_SIZE($input));
   else {
    PyErr_Clear();
    PyErr_SetString(PyExc_TypeError,"not a bytes or bytearray");
    return NULL;
   }
}

%typemap(freearg) const std::string& {
  delete $1;
}

namespace GeoCal {
std::string serialize_write_binary(const boost::shared_ptr<GenericObject>& Obj);

boost::shared_ptr<GenericObject> 
serialize_read_binary(const std::string& Data);
}

// List of things "import *" will include
%python_export("have_serialize_supported","serialize_write","serialize_write_string","serialize_read_generic","serialize_read_generic_string","serialize_write_binary","serialize_read_binary", "serialize_read_binary_generic")
