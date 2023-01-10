%include "fp_common.i"

%{
#include "multiscatt_rt_driver.h"
%}

%base_import(spurr_rt_driver)

%import "spurr_interface_types.i"
%import "spurr_interface_masters.i"

%fp_shared_ptr(FullPhysics::MultiScattRtDriver);

namespace FullPhysics {

class MultiScattRtDriver : public SpurrRtDriver {
public:
  MultiScattRtDriver(bool do_solar = true, bool do_thermal = false);

  %python_attribute(number_moment, int)
  %python_attribute(number_stream, int)
  %python_attribute(surface_type, int)
  %python_attribute(do_multi_scatt_only, bool)
  %python_attribute(do_thermal_scattering, bool)
  %python_attribute(pure_nadir, bool)

  virtual void notify_update(const RtAtmosphere& atm);

  void setup_sphericity(blitz::Array<double, 1> zen, bool do_multi_scatt_only, bool pure_nadir);
  void set_plane_parallel();
  void set_pseudo_spherical();
  void set_plane_parallel_plus_ss_correction();
  void set_line_of_sight();

  void setup_height_grid(const blitz::Array<double, 1>& height_grid);
  void setup_geometry(double sza, double azm, double zen);

  void setup_thermal_inputs(double surface_bb, const blitz::Array<double, 1> atmosphere_bb);

  void setup_optical_inputs(const blitz::Array<double, 1>& od, 
                            const blitz::Array<double, 1>& ssa,
                            const blitz::Array<double, 3>& pf);

  void clear_linear_inputs();
  void setup_linear_inputs(const ArrayAd<double, 1>& od,
                           const ArrayAd<double, 1>& ssa,
                           const ArrayAd<double, 3>& pf,
                           bool do_surface_linearization);
  
  void calculate_rt() const;

  double get_intensity() const;
  void copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf_param, double& jac_surf_temp, blitz::Array<double, 1>& jac_atm_temp) const;

  %python_attribute_abstract(rt_interface, boost::shared_ptr<Spurr_Lps_Masters_Base>)
  %python_attribute_abstract(brdf_interface, boost::shared_ptr<Spurr_Brdf_Lin_Sup_Masters_Base>)

  %pickle_serialization();
};

}
