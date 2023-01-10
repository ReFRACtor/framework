%include "fp_common.i"

%{
#include "spurr_rt_driver.h"
%}

%base_import(generic_object)

%import "spurr_brdf_driver.i"

%import "array_ad.i"
%import "rt_atmosphere.i"

%fp_shared_ptr(FullPhysics::SpurrRtDriver);

namespace FullPhysics {

class SpurrRtDriver : public GenericObject {
public:

  SpurrRtDriver(bool do_solar = true, bool do_thermal = false);

  std::string print_to_string() const;
  virtual void notify_update(const RtAtmosphere& atm);
  virtual double reflectance_calculate(const blitz::Array<double, 1>& height_grid,
                                       double sza, double azm, double zen,
                                       int surface_type,
                                       const blitz::Array<double, 1>& surface_parameters,
                                       const blitz::Array<double, 1>& od, 
                                       const blitz::Array<double, 1>& ssa,
                                       const blitz::Array<double, 3>& pf,
                                       double surface_bb = 0,
                                       const blitz::Array<double, 1>& atmosphere_bb = blitz::Array<double,1>());

  virtual void reflectance_and_jacobian_calculate(const blitz::Array<double, 1>& height_grid,
                                                  double sza, double azm, double zen,
                                                  int surface_type,
                                                  ArrayAd<double, 1>& surface_parameters,
                                                  const ArrayAd<double, 1>& od, 
                                                  const ArrayAd<double, 1>& ssa,
                                                  const ArrayAd<double, 3>& pf,
                                                  double& reflectance,
                                                  blitz::Array<double, 2>& jac_atm, 
                                                  blitz::Array<double, 1>& jac_surf_param,
                                                  double &jac_surf_temp,
                                                  blitz::Array<double, 1>& jac_atm_temp,
                                                  double surface_bb = 0,
                                                  const blitz::Array<double, 1>& atmosphere_bb = blitz::Array<double,1>());

  %python_attribute(brdf_driver, boost::shared_ptr<SpurrBrdfDriver>)
  virtual void setup_height_grid(const blitz::Array<double, 1>& height_grid)  = 0;
  virtual void setup_geometry(double sza, double azm, double zen) = 0;
  virtual void setup_thermal_inputs(double surface_bb, const blitz::Array<double, 1> atmosphere_bb) = 0;
  virtual void setup_optical_inputs(const blitz::Array<double, 1>& od, 
                                    const blitz::Array<double, 1>& ssa,
                                    const blitz::Array<double, 3>& pf) = 0;
  virtual void clear_linear_inputs()  =  0;
  virtual void setup_linear_inputs(const ArrayAd<double, 1>& od,
                                   const ArrayAd<double, 1>& ssa,
                                   const ArrayAd<double, 3>& pf,
                                   bool do_surface_linearization) = 0;
  virtual void calculate_rt() const = 0;
  virtual double get_intensity() const = 0;
  virtual void copy_jacobians(blitz::Array<double, 2>& jac_atm, blitz::Array<double, 1>& jac_surf_params, double& jac_surf_temp, blitz::Array<double, 1>& jac_atm_temp) const = 0;
  %pickle_serialization();
};

}
