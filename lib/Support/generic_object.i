// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"
%{
#include "generic_object.h"
#include <cstdint>  
%}

%shared_ptr(FullPhysics::GenericObject)

%typemap(out) const boost::shared_ptr< FullPhysics::GenericObject >& {
  %set_output(FullPhysics::swig_to_python_or_none($1));
}

%typemap(out) boost::shared_ptr< FullPhysics::GenericObject >& {
  %set_output(FullPhysics::swig_to_python_or_none($1));
}

%typemap(out) boost::shared_ptr< FullPhysics::GenericObject > {
  %set_output(FullPhysics::swig_to_python_or_none($1));
}

namespace FullPhysics {
class GenericObject {
public:
  virtual ~GenericObject();
  %extend {
//-----------------------------------------------------------------------
/// We generally convert object returned to the most specific object
/// in python. However we may have instances where this conversion
/// doesn't already happen (e.g., we have a
/// std::vector<boost::shared<MoreGeneral> >). This function uses the
/// SWIG magic to convert this if needed.
//-----------------------------------------------------------------------
    static boost::shared_ptr<GenericObject> convert_to_most_specific_class
      (const boost::shared_ptr<GenericObject>& V) { return V; }

//-----------------------------------------------------------------------
/// Because the swig wrapper is different each time we get an object,
/// it can be hard to tell if two objects are the same. This returns
/// true if V points to the same object.    
//-----------------------------------------------------------------------
    
    bool is_same_ptr(const boost::shared_ptr<GenericObject>& V) const
    { return $self == V.get(); }
  }
};
}

%template(vector_GenericObject) std::vector<boost::shared_ptr<FullPhysics::GenericObject> >;

// List of things "import *" will include
%python_export("GenericObject", "vector_GenericObject")
