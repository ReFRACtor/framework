%include "fp_common.i"

%{
#include "optical_properties_imp_base.h"
#include "sub_state_vector_array.h"
%}

%base_import(optical_properties)

%import "array_ad.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::AerosolPhaseFunctionHelper)
%fp_shared_ptr(FullPhysics::AerosolPhaseFunctionPassThruHelper)
%fp_shared_ptr(FullPhysics::AerosolPhaseFunctionComputeHelper)

%fp_shared_ptr(FullPhysics::OpticalPropertiesImpBase)

%include "optical_properties_imp_base.h"
