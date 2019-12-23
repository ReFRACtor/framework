%include "fp_common.i"

%{
#include "optical_properties_lsi.h"
#include "sub_state_vector_array.h"
%}

%base_import(optical_properties_imp_base)

%import "aerosol_optical.i"

%fp_shared_ptr(FullPhysics::OpticalPropertiesLsi)

%include "optical_properties_lsi.h"
