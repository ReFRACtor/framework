%include "fp_common.i"
%{
#include "pca_binning.h"
#include "sub_state_vector_array.h"
#include "altitude.h"
%}

%import "optical_properties.i"

%fp_shared_ptr(FullPhysics::PCABinning);

%include "pca_binning.h"
