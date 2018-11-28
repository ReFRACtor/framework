%include "common.i"
%{
#include "pca_optical_properties.h"
#include "sub_state_vector_array.h"
%}

%import "atmosphere_oco.i"
%import "spectral_domain.i"
%import "hdf_file.i"

%fp_shared_ptr(FullPhysics::PCAOpticalProperties);
%fp_shared_ptr(FullPhysics::PCAOpticalPropertiesAtmosphere);
%fp_shared_ptr(FullPhysics::PCAOpticalPropertiesFile);

%include "pca_optical_properties.h"
