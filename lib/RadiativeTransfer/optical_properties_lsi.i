// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"

%{
#include "optical_properties_lsi.h"
#include "sub_state_vector_array.h"
#include "altitude.h"
%}

%base_import(optical_properties_imp_base)

%import "aerosol_optical.i"
%import "optical_properties_wrt_rt.i"

%fp_shared_ptr(FullPhysics::OpticalPropertiesLsi)

%template (vector_optical_properties_lsi) std::vector<boost::shared_ptr<FullPhysics::OpticalPropertiesLsi> >;

namespace FullPhysics {
class OpticalPropertiesLsi : public OpticalPropertiesImpBase {
public:
  OpticalPropertiesLsi(const ArrayAd<double, 2>& packed_properties, double wavenumber, const boost::shared_ptr<AerosolOptical>& aerosol, int num_gas, int num_aerosol);
  static ArrayAd<double, 2> pack(const boost::shared_ptr<OpticalPropertiesWrtRt>& source_properties);
  %pickle_serialization();
};

}

