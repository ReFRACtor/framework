// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// This wraps a shared_ptr up for FullPhysics. In particular, we use
// RTTI to make sure that the most specific type is returned in
// python, rather than the most general type. This maps better to the
// standard duck typing done in python. See swig_cast_test.py for a
// simple example of this.

%{
#define SWIG_MAPPER_NAMESPACE FullPhysics
%}

%include "shared_ptr_type_mapper.i"

%include "shared_ptr_type_mapper.i"

%define %fp_shared_ptr(TYPE...)
%shared_ptr_type_mapper(FullPhysics, TYPE)
%enddef

