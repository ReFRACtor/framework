// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

//--------------------------------------------------------------
// This provides common functions etc. used throughout our SWIG
// code. This should get included at the top of each swig file.
//--------------------------------------------------------------

// The module actually gets overridden by SWIG command line options
// when we build. But we still need to supply this to get the
// directors=1 and allprotected=1 set.

%module(directors="1", allprotected="1") geocal

%{
#include <boost/shared_ptr.hpp>
%}

// Short cut for ingesting a base class
%define %base_import(NAME)
%import(module="skeleton_swig.NAME") "NAME.i"
%enddef

// Map std::string to and from the native string type
%naturalvar std::string;

%include <std_vector.i>

// Include our own rules and common imports here.
%include "skeleton_shared_ptr.i"
%include "swig_exception.i"
%include "swig_print.i"
%include "swig_python_attribute.i"
%include "swig_pickle.i"
%import "swig_std.i"
%include "swig_array_inc.i"
%import "swig_array.i"
%include "swig_boost_array_inc.i"
%import "swig_boost_array.i"
%import "swig_boost_optional.i"
%include "swig_iostream_inc.i"
%import "swig_iostream.i"

%{

//--------------------------------------------------------------
/// The default conversion of a python sequence to a std::vector<T>
/// doesn't work correctly for boost::shared_ptr types. This
/// expression here does this correctly.
//--------------------------------------------------------------

namespace swig {
  template <class T>
  struct traits_asptr<std::vector<boost::shared_ptr<T> > >  {
    static int asptr(PyObject *obj, 
		     std::vector<boost::shared_ptr<T> > **vec) {
      // Long name, so shorten
      typedef std::vector<boost::shared_ptr<T> > vtype;
      if (obj == Py_None || SWIG_Python_GetSwigThis(obj)) {
	vtype *p;
	if (::SWIG_ConvertPtr(obj,(void**)&p,
			      swig::type_info<vtype>(),0) == SWIG_OK) {
	  if (vec) *vec = p;
	  return SWIG_OLDOBJ;
	}
      } else if (PySequence_Check(obj)) {
	try {
	  if (vec) {
	    vtype *pseq = new vtype();
	    PyObject *iterator = PyObject_GetIter(obj);
	    while(PyObject *item = PyIter_Next(iterator)) {
	      boost::shared_ptr<T> *itemp;
	      int newmem = 0;
	      int res = SWIG_ConvertPtrAndOwn(item, (void**) &itemp, 
			swig::type_info<boost::shared_ptr<T> >(),  0, &newmem);
	      if(!SWIG_IsOK(res)) {
		Py_DECREF(item);
		Py_DECREF(iterator);
		return SWIG_ERROR;
	      }
	      pseq->push_back(*itemp);
	      Py_DECREF(item);
	    }
	    Py_DECREF(iterator);
	    *vec = pseq;
	    return SWIG_NEWOBJ;
	  } else {
	    SwigPySequence_Cont<boost::shared_ptr<T> > swigpyseq(obj);
	    return swigpyseq.check() ? SWIG_OK : SWIG_ERROR;
	  }
	} catch (std::exception& e) {
	  if (vec) {
	    if (!PyErr_Occurred()) {
	      PyErr_SetString(PyExc_TypeError, e.what());
	    }
	  }
	  return SWIG_ERROR;
	}
      }
      return SWIG_ERROR;
    }
  };
}
%}
