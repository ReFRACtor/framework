%include "common.i"
%{
#include "pca_optical_properties.h"
#include "sub_state_vector_array.h"
%}

%import "atmosphere_oco.i"
%import "spectral_domain.i"

%fp_shared_ptr(FullPhysics::PCAOpticalProperties);
%fp_shared_ptr(FullPhysics::PCAOpticalPropertiesAtmosphere);

%include "pca_optical_properties.h"
