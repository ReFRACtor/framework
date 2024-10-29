// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "fp_exception.h"
%}

%base_import(generic_object)

%init {
  FullPhysics::no_gsl_abort();
}

// Rename to avoid conflict with built in Python Exception object
%rename(FpException) Exception;

%fp_shared_ptr(FullPhysics::Exception)
namespace FullPhysics {
  class Exception : public GenericObject {
  public:
    Exception(const std::string& W);
    std::string print_to_string() const;
    std::string print_parent() const;
    const char* what();
  };

void no_gsl_abort();
}

