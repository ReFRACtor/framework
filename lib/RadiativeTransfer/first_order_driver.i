// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"
%{
#include "first_order_driver.h"
%}

%import "first_order_interface.i"

%base_import(spurr_driver)

%fp_shared_ptr(FullPhysics::FirstOrderDriver);

namespace FullPhysics {

class FirstOrderDriver : public SpurrRtDriver {
public:
  FirstOrderDriver(int number_layers, int surface_type, int number_streams, int number_moments,
		   bool do_solar = true, bool do_thermal = false);
  %python_attribute(number_moment, int);
  %python_attribute(number_stream, int);
  %python_attribute(number_layer, int);
  %python_attribute(surface_type, int);
  void set_plane_parallel();
  void set_pseudo_spherical();
  void setup_height_grid(const blitz::Array<double, 1>& height_grid) const;
  void setup_geometry(double sza, double azm, double zen) const;
  void setup_thermal_inputs(double surface_bb,
			    const blitz::Array<double, 1>& atmosphere_bb) const;
  void setup_optical_inputs(const blitz::Array<double, 1>& od, 
			    const blitz::Array<double, 1>& ssa,
			    const blitz::Array<double, 2>& pf) const;
  void clear_linear_inputs() const;
  void setup_linear_inputs(const ArrayAd<double, 1>& od,
			   const ArrayAd<double, 1>& ssa,
			   const ArrayAd<double, 2>& pf,
			   bool do_surface_linearization) const;
  void calculate_rt() const;
  double get_intensity() const;
  void copy_jacobians(blitz::Array<double, 2>& jac_atm,
		      blitz::Array<double, 1>& jac_surf_param,
		      double& jac_surf_temp,
		      blitz::Array<double, 1>& jac_atm_temp) const;
  const boost::shared_ptr<Fo_Ssgeometry_Master> geometry_interface() const;
  const boost::shared_ptr<Fo_Scalarss_Rtcalcs_Ilps> solar_interface() const;
  %pickle_serialization();
};
}

