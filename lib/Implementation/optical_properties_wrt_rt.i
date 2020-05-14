%include "fp_common.i"

%{
#include "optical_properties_wrt_rt.h"
#include "sub_state_vector_array.h"
#include "altitude.h"
%}

%base_import(optical_properties_init_base)

%fp_shared_ptr(FullPhysics::OpticalPropertiesWrtRt)

%template (vector_optical_properties_wrt_rt) std::vector<boost::shared_ptr<FullPhysics::OpticalPropertiesWrtRt> >;

%include "optical_properties_wrt_rt.h"
