%include "common.i"
%{
#include "pca_eigensolver.h"
%}

%fp_shared_ptr(FullPhysics::PCAEigenSolver);

%include "pca_eigensolver.h"
