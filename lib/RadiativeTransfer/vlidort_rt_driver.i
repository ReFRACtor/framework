%include "fp_common.i"

%{
#include "vlidort_rt_driver.h"
%}

%base_import(multiscatt_rt_driver)

%import "vlidort_interface_masters.i"
%import "vlidort_interface_types.i"
%import "vlidort_brdf_driver.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::VLidortRtDriver);

namespace FullPhysics {

class VLidortRtDriver : public MultiScattRtDriver {
public:
  VLidortRtDriver(int nstream, int nmoment, bool do_multi_scatt_only, int surface_type,
          const blitz::Array<double, 1>& zen, bool pure_nadir, int nstokes = 3,
          bool do_solar_sources = true, bool do_thermal_emission = false, bool do_thermal_scattering = true);

  %python_attribute(number_stokes, int)

  %python_attribute(rt_interface, boost::shared_ptr<Spurr_Lps_Masters_Base>)
  %python_attribute(brdf_interface, boost::shared_ptr<Spurr_Brdf_Lin_Sup_Masters_Base>)

  %python_attribute(vlidort_brdf_driver, boost::shared_ptr<VLidortBrdfDriver>)
  %python_attribute(vlidort_brdf_interface, boost::shared_ptr<VBrdf_Linsup_Masters>)
  %python_attribute(vlidort_interface, boost::shared_ptr<VLidort_Lps_Masters>)

  const blitz::Array<double, 1> get_intensity() const;
  void copy_jacobians(blitz::Array<double, 3>& jac_atm, blitz::Array<double, 2>& jac_surf_param, blitz::Array<double, 1>& jac_surf_temp, blitz::Array<double, 2>& jac_atm_temp) const;

  %pickle_serialization();
};

}
