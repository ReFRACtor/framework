// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// This wraps a shared_ptr up for a given namespace. In particular, we use
// RTTI to make sure that the most specific type is returned in
// python, rather than the most general type. This maps better to the
// standard duck typing done in python.

// There appears to be a bug in the shared_ptr handler of SWIG, as of
// version 3.0.12. This is the normal handler, with some fixes added.
%include "my_shared_ptr.i"

%{
#include "swig_type_mapper.h"
#include <boost/make_shared.hpp>
%}

%define %shared_ptr_type_mapper(NAMESPACE, TYPE...)
%shared_ptr(TYPE)
%init {
  NAMESPACE::SwigTypeMapperBase::add(typeid(TYPE), boost::make_shared<NAMESPACE::SwigTypeMapper< TYPE > > ("boost::shared_ptr< TYPE > *"));
}

%typemap(out) const boost::shared_ptr< TYPE >& {
  %set_output(NAMESPACE::swig_to_python($1));
}

%typemap(out) boost::shared_ptr< TYPE >& {
  %set_output(NAMESPACE::swig_to_python($1));
}

%typemap(out) boost::shared_ptr< TYPE > {
  %set_output(NAMESPACE::swig_to_python($1));
}

%typemap(out) boost::shared_ptr< TYPE >* {
  %set_output(NAMESPACE::swig_to_python(*$1));
}

%typemap(out) const boost::shared_ptr< TYPE > * {
  %set_output(NAMESPACE::swig_to_python(*$1));
}

%typemap(out) const boost::shared_ptr< TYPE > *& {
  %set_output(NAMESPACE::swig_to_python(*$1));
}

%enddef

