%include "common.i"
%{
#include "pca_binning.h"
#include "sub_state_vector_array.h"
%}

%import "atmosphere_oco.i"
%import "spectral_domain.i"

%fp_shared_ptr(FullPhysics::PCABinning);

%include "pca_binning.h"
