%include "fp_common.i"

%{
#include "optical_properties.h"
#include "sub_state_vector_array.h"
%}

%base_import(generic_object)

%import "array_ad.i"
%import "double_with_unit.i"
%import "absorber.i"
%import "rayleigh.i"
%import "aerosol.i"

%fp_shared_ptr(FullPhysics::OpticalPropertiesImpBase)
%fp_shared_ptr(FullPhysics::OpticalPropertiesWrtInput)
%fp_shared_ptr(FullPhysics::OpticalPropertiesWrtRt)

%include "optical_properties.h"
