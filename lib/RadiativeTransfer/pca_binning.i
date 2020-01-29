%include "fp_common.i"
%{
#include "pca_binning.h"
#include "sub_state_vector_array.h"
%}

%import "optical_properties_wrt_rt.i"

%fp_shared_ptr(FullPhysics::PCABinning);

%include "pca_binning.h"
