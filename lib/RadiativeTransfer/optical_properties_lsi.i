%include "fp_common.i"

%{
#include "optical_properties_lsi.h"
#include "sub_state_vector_array.h"
#include "altitude.h"
%}

%base_import(optical_properties_imp_base)

%import "aerosol_optical.i"
%import "optical_properties_wrt_rt.i"

%fp_shared_ptr(FullPhysics::OpticalPropertiesLsi)

%template (vector_optical_properties_lsi) std::vector<boost::shared_ptr<FullPhysics::OpticalPropertiesLsi> >;

%include "optical_properties_lsi.h"
