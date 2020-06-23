// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

//--------------------------------------------------------------
// Support for writing the __reduce__ function needed for pickling.
//--------------------------------------------------------------

%{
#include "serialize_function.h"
#include <stdexcept>
// This is defined in swig_wrap.tmpl, so it gets put into swig_wrap.cc
std::string parse_python_exception();
%}

%{
//--------------------------------------------------------------
/// Support routines for calling cPickle.dumps from C
//--------------------------------------------------------------

inline PyObject* cpickle_module()
{
  static PyObject* mod = 0;
  if(!mod)
    mod = PyImport_ImportModule("pickle");
  return mod;
}

inline std::string cpickle_dumps(PyObject* obj)
{
  PyObject* res = PyObject_CallMethodObjArgs(cpickle_module(),
					     PyString_FromString("dumps"),
					     obj, NULL);
  if(PyErr_Occurred()) {
    throw std::runtime_error("Python error occurred:\n" + parse_python_exception());
  }
  return std::string(PyString_AsString(res));
}

inline PyObject* cpickle_loads(const std::string& S)
{
  PyObject* res = PyObject_CallMethodObjArgs(cpickle_module(),
					     PyString_FromString("loads"),
					     PyString_FromString(S.c_str()), 
					     NULL);
  if(PyErr_Occurred()) {
    throw std::runtime_error("Python error occurred:\n" + parse_python_exception());
  }
  return res;
}

%}
//--------------------------------------------------------------
// Code to support the python side.
//--------------------------------------------------------------

%pythoncode {
import os

def _new_from_init(cls, version, *args):
    '''For use with pickle, covers common case where we just store the
    arguments needed to create an object. See for example HdfFile'''
    if(cls.pickle_format_version() != version):
      raise RuntimeException("Class is expecting a pickled object with version number %d, but we found %d" % (cls.pickle_format_version(), version))
    inst = cls.__new__(cls)
    inst.__init__(*args)
    return inst
 
def _new_from_serialization(data):
    return SWIG_MODULE.serialize_read_binary(data)

def _new_from_serialization_dir(dir, data):
    curdir = os.getcwd()
    try:
      os.chdir(dir)
      return SWIG_MODULE.serialize_read_binary(data)
    finally:
      os.chdir(curdir)
	
	
def _new_vector(cls, version, lst):
    '''Create a vector from a list.'''
    if(cls.pickle_format_version() != version):
      raise RuntimeException("Class is expecting a pickled object with version number %d, but we found %d" % (cls.pickle_format_version(), version))
    inst = cls.__new__(cls)
    inst.__init__()
    for i in lst:
       inst.append(i)
    return inst

def _new_from_set(cls, version, *args):
    '''For use with pickle, covers common case where we use a set function 
    to assign the value'''
    if(cls.pickle_format_version() != version):
      raise RuntimeException("Class is expecting a pickled object with version number %d, but we found %d" % (cls.pickle_format_version(), version))
    inst = cls.__new__(cls)
    inst.__init__()
    inst.set(*args)
    return inst
}

%define %pickle_init(VER, ARG...)
  %pythoncode {
@classmethod
def pickle_format_version(cls):
  return VER

def __reduce__(self):
  return _new_from_init, (self.__class__, VER, ARG)
}
%enddef

%define %pickle_serialization()
  %pythoncode {
def __reduce__(self):
  return _new_from_serialization, (SWIG_MODULE.serialize_write_binary(self),)
}
%enddef

%define %pickle_serialization_dir()
  %pythoncode {
def __reduce__(self):
  return _new_from_serialization_dir, (os.getcwd(), SWIG_MODULE.serialize_write_binary(self),)
}
%enddef

%define %pickle_vector()
  %pythoncode {
@classmethod
def pickle_format_version(cls):
  return 1

def to_list(self):
   res = []
   for i in range(self.size()):
      res.append(self[i])
   return res

def __reduce__(self):
  return _new_vector, (self.__class__, 1, self.to_list())
}
%enddef
