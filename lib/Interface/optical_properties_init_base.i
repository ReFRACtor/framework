// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "optical_properties_init_base.h"
#include "sub_state_vector_array.h"
#include "altitude.h"
%}

%base_import(optical_properties_imp_base)

%import "absorber.i"
%import "rayleigh.i"
%import "aerosol.i"

%fp_shared_ptr(FullPhysics::OpticalPropertiesInitBase);
namespace FullPhysics {
class OpticalPropertiesInitBase : public OpticalPropertiesImpBase {
public:
  OpticalPropertiesInitBase();
  virtual void initialize(const ArrayAd<double, 1>& rayleigh_od, 
                          const ArrayAd<double, 2>& gas_od,
                          const int num_jacobians = -1);
  virtual void initialize(const ArrayAd<double, 1>& rayleigh_od, 
                          const ArrayAd<double, 2>& gas_od,
                          const ArrayAd<double, 2>& aerosol_ext_od,
                          const ArrayAd<double, 2>& aerosol_sca_od,
                          const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments,
                          const int num_jacobians = -1);
  virtual void initialize(const DoubleWithUnit spectral_point,
                          const int channel_index,
                          const boost::shared_ptr<Absorber>& absorber,
                          const boost::shared_ptr<Rayleigh>& rayleigh,
                          const boost::shared_ptr<Aerosol>& aerosol,
                          const int num_jacobians = -1);
  %pickle_serialization();
protected:
  virtual void initialize_with_jacobians(const ArrayAd<double, 1>& rayleigh_od, 
                                         const ArrayAd<double, 2>& gas_od,
                                         const ArrayAd<double, 2>& aerosol_ext_od,
                                         const ArrayAd<double, 2>& aerosol_sca_od,
                                         const boost::shared_ptr<AerosolPhaseFunctionHelper>& aer_pf_helper,
                                         const int num_jacobians = -1) = 0;
};

}

