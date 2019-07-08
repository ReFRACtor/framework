%include "fp_common.i"

%{
#include "jacobian_size_update.h"
%}

%base_import(state_vector_observer)

%import "jacobian_size_mixin.i"

%fp_shared_ptr(FullPhysics::JacobianSizeUpdate);

%include "jacobian_size_update.h"
