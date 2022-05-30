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
%include "swig_vector_shared_ptr.i"
