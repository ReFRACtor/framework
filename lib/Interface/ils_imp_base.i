%include "common.i"

%{
#include "ils_imp_base.h"
#include "sub_state_vector_array.h"
%}

%import "double_with_unit.i"
%import "array_ad.i"
%import "spectral_domain.i"
%import "double_with_unit.i"

%base_import(ils)
%base_import(observer)
%base_import(sample_grid)

%fp_shared_ptr(FullPhysics::IlsImpBase);

%feature("director") FullPhysics::IlsImpBase;

%include "ils_imp_base.h"
