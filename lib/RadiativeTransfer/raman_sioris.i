%include "fp_common.i"

%{
#include "raman_sioris.h"

// These are needed because of the inclusion of AtmosphereStandard
#include "instrument.h"
#include "spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
%}

%base_import(pressure)
%base_import(spectrum_effect_imp_base)

%import "spectral_domain.i"
%import "double_with_unit.i"
%import "atmosphere_standard.i"
%import "solar_model.i"

%fp_shared_ptr(FullPhysics::RamanSiorisEffect);

%include "raman_sioris.h"
