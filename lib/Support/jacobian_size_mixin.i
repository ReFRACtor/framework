%include "fp_common.i"

%{
#include "jacobian_size_mixin.h"
%}

%fp_shared_ptr(FullPhysics::JacobianSizeMixin);

%include "jacobian_size_mixin.h"
