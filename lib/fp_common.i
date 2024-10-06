// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

//--------------------------------------------------------------
// This provides common functions etc. used throughout our SWIG
// code. This should get included at the top of each swig file.
//--------------------------------------------------------------

// The module actually gets overridden by SWIG command line options
// when we build. But we still need to supply this to get the
// directors=1 and allprotected=1 set.
//
// See https://github.com/swig/swig/issues/2260 for discussion of
// moduleimport, this is needed because we place everything in one
// library _swig_wrap.so. Swig 4 seems to prefer splitting this into
// separate libraries (so for example "exception.py" would import
// _exception.so). It is possible we can change to that in the future,
// although it seems like it is useful to have it all one library like
// things like cython do. But for now we can just work around this.
// The effect of this is just to change the "import" line in the swig
// generated python code.

#if SWIG_VERSION < 0x040000  
%module(directors="1", allprotected="1") foo
#else
%module(moduleimport="from ._swig_wrap import $module", directors="1", allprotected="1") foo
#endif
#define SWIG_MODULE_ALREADY_DONE 1
#define SWIG_MODULE refractor.framework_swig

// Starting with swig 4.2, we get messages about memory leaks of
// pointers to some types. Not really clear what is going on
// here. but we really just want this silenced
%begin %{
// Don't report memory leaks that we can't do anything about
#define SWIG_PYTHON_SILENT_MEMLEAK 1  
%}

%{
#include <boost/shared_ptr.hpp>
#include <boost/rational.hpp>
%}

// Short cut for ingesting a base class
%define %base_import(NAME)
%import(module="refractor.framework_swig.NAME") "NAME.i"
%enddef

// Map std::string to and from the native string type
%naturalvar std::string;

%include <std_vector.i>

// Include our own rules and common imports here.
%include "fp_shared_ptr.i"
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
%include "swig_vector_shared_ptr.i"
%import "swig_rational.i"
%include "director_serialization.i"
%define %fp_director_serialization(BNAME, TYPE...)
%director_serialization(FullPhysics, fp_serialize_support.h, BNAME, TYPE)
%enddef
