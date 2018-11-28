%include "common.i"
%{
#include "pca_eigensolver.h"
%}

%fp_shared_ptr(FullPhysics::PCAEigenSolver);
%fp_shared_ptr(FullPhysics::PCAEigenSolverGeneric);
%fp_shared_ptr(FullPhysics::PCAEigenSolverFortran);

%include "pca_eigensolver.h"
