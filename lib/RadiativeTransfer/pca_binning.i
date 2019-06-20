%include "fp_common.i"
%{
#include "pca_binning.h"
%}

%import "pca_optical_properties.i"

%fp_shared_ptr(FullPhysics::PCABinning);

%include "pca_binning.h"
