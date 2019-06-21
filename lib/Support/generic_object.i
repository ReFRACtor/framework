// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"
%{
#include "generic_object.h"
%}

%shared_ptr(FullPhysics::GenericObject)
namespace FullPhysics {
class GenericObject {
public:
  virtual ~GenericObject();
};
}

%typemap(out) const boost::shared_ptr< FullPhysics::GenericObject >& {
  %set_output(FullPhysics::swig_to_python_or_none($1));
}

%typemap(out) boost::shared_ptr< FullPhysics::GenericObject >& {
  %set_output(FullPhysics::swig_to_python_or_none($1));
}

%typemap(out) boost::shared_ptr< FullPhysics::GenericObject > {
  %set_output(FullPhysics::swig_to_python_or_none($1));
}

// List of things "import *" will include
%python_export("GenericObject")
