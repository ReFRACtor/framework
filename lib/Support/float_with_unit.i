// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "float_with_unit.h"
%}

%base_import(generic_object)
%import "unit.i"
%fp_shared_ptr(FullPhysics::FloatWithUnit)
namespace FullPhysics {
using FloatWithUnit = ScalarWithUnit<float>;
}

