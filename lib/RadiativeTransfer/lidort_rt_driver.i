%include "fp_common.i"

%{
#include "lidort_rt_driver.h"
%}

%base_import(multiscatt_rt_driver)

%import "lidort_interface_masters.i"
%import "lidort_interface_types.i"
%import "lidort_brdf_driver.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::LidortRtDriver);

namespace FullPhysics {

class LidortRtDriver : public MultiScattRtDriver {
public:
  LidortRtDriver(int nstream, int nmoment, bool do_multi_scatt_only, int surface_type, 
          const blitz::Array<double, 1>& zen, bool pure_nadir, 
          bool do_solar_sources = true, bool do_thermal_emission = false, bool do_thermal_scattering = true);

  void set_plane_parallel();
  void set_pseudo_spherical();
  void set_plane_parallel_plus_ss_correction();
  void set_line_of_sight();

  %python_attribute(rt_interface, boost::shared_ptr<Spurr_Lps_Masters_Base>)
  %python_attribute(brdf_interface, boost::shared_ptr<Spurr_Brdf_Lin_Sup_Masters_Base>)

  %python_attribute(lidort_brdf_driver, boost::shared_ptr<LidortBrdfDriver>)
  %python_attribute(lidort_brdf_interface, boost::shared_ptr<Brdf_Lin_Sup_Masters>)
  %python_attribute(lidort_interface, boost::shared_ptr<Lidort_Lps_Masters>)

  double get_intensity() const;
  void copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf_param, double& jac_surf_temp, blitz::Array<double, 1>& jac_atm_temp) const;

  %pickle_serialization();
};

}
