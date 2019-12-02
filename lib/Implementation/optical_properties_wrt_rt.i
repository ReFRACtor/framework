%include "fp_common.i"

%{
#include "optical_properties_wrt_rt.h"
#include "sub_state_vector_array.h"
%}

%base_import(optical_properties_imp_base)

%fp_shared_ptr(FullPhysics::OpticalPropertiesWrtRt)

%include "optical_properties_wrt_rt.h"
