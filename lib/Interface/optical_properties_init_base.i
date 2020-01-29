%include "fp_common.i"

%{
#include "optical_properties_init_base.h"
#include "sub_state_vector_array.h"
%}

%base_import(optical_properties_imp_base)

%import "absorber.i"
%import "rayleigh.i"
%import "aerosol.i"

%fp_shared_ptr(FullPhysics::OpticalPropertiesInitBase)

%include "optical_properties_init_base.h"
