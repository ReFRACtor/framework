%include "fp_common.i"

%{
#include "optical_properties_pca.h"
#include "sub_state_vector_array.h"
%}

%base_import(optical_properties_imp_base)

%import "aerosol_optical.i"
%import "optical_properties_wrt_rt.i"

%fp_shared_ptr(FullPhysics::OpticalPropertiesPca)

%template (vector_optical_properties_pca) std::vector<boost::shared_ptr<FullPhysics::OpticalPropertiesPca> >;

%include "optical_properties_pca.h"
