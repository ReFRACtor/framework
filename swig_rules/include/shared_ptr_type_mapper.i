// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// This wraps a shared_ptr up for a given namespace. In particular, we use
// RTTI to make sure that the most specific type is returned in
// python, rather than the most general type. This maps better to the
// standard duck typing done in python.

// The SWIG handler for shared_ptr doesn't properly handle directors
// that are python classes (see swig_rules/lib/DirectorNotes.md. This
// is our version of share_ptr.i that adds proper handling of this.
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

