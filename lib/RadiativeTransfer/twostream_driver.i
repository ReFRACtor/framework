// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "twostream_driver.h"
%}
%base_import(spurr_driver)
%import "twostream_interface.i"
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::TwostreamBrdfDriver);
%fp_shared_ptr(FullPhysics::TwostreamRtDriver);

namespace FullPhysics {
class TwostreamBrdfDriver : public SpurrBrdfDriver {
public:
  TwostreamBrdfDriver(int surface_type);
  virtual ~TwostreamBrdfDriver();

  virtual void setup_geometry(double sza, double azm, double zen);

  %python_attribute(n_brdf_kernels,int)
  %python_attribute(n_kernel_factor_wfs,int)
  %python_attribute(n_kernel_params_wfs,int)
  %python_attribute(n_surface_wfs,int)
  %python_attribute(do_shadow_effect,bool)
  %python_attribute(brdf_interface, boost::shared_ptr<Twostream_Ls_Brdf_Supplement>)
  virtual bool do_kparams_derivs(const int kernel_index) const;
  %pickle_serialization();
};

class TwostreamRtDriver : public SpurrRtDriver {
public:
  TwostreamRtDriver(int nlayers, int surface_type, bool do_fullquadrature = true,
          bool do_solar = true, bool do_thermal = false);

  void setup_height_grid(const blitz::Array<double, 1>& height_grid);
  void setup_geometry(double sza, double azm, double zen);
  void setup_thermal_inputs(double surface_bb, const blitz::Array<double, 1> atmosphere_bb);
  void setup_optical_inputs(const blitz::Array<double, 1>& od, 
                            const blitz::Array<double, 1>& ssa,
                            const blitz::Array<double, 2>& pf);
  void clear_linear_inputs();
  void setup_linear_inputs(const ArrayAd<double, 1>& od,
                           const ArrayAd<double, 1>& ssa,
                           const ArrayAd<double, 2>& pf,
                           bool do_surface_linearization);

  void calculate_rt() const;
  double get_intensity() const;
  void copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf_param, double& jac_surf_temp, blitz::Array<double, 1>& jac_atm_temp) const;

  %python_attribute(twostream_brdf_driver, boost::shared_ptr<TwostreamBrdfDriver>)
  %python_attribute(brdf_interface, boost::shared_ptr<Twostream_Ls_Brdf_Supplement>)
  %python_attribute(twostream_interface, boost::shared_ptr<Twostream_Lps_Master>)

  %python_attribute(do_full_quadrature, bool)
  %pickle_serialization();
};
}
