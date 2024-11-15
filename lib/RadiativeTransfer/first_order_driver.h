#ifndef FO_DRIVER_H
#define FO_DRIVER_H

#include "spurr_driver.h"
#include "first_order_interface.h"

// Include BRDF driver interface from LIDORT driver
#include "lidort_driver.h"

/****************************************************************//**
  Driver for optimized first order of scattering RT
 *******************************************************************/

namespace FullPhysics {

class FirstOrderDriver : public SpurrRtDriver {
public:
  FirstOrderDriver(int number_layers, int surface_type, int number_streams, int number_moments,
                   bool do_solar = true, bool do_thermal = false); 

  virtual void notify_update(const RtAtmosphere& atm);

  int number_moment() const { return num_moments_; }
  int number_stream() const { return num_streams_; }
  int number_layers() const { return num_layers_; }
  int surface_type() const { return surface_type_; }

  void set_plane_parallel();
  void set_pseudo_spherical();
  void set_line_of_sight();
  
  void setup_height_grid(const blitz::Array<double, 1>& height_grid);
  void setup_geometry(double sza, double azm, double zen);
  
  void setup_thermal_inputs(double surface_bb, const blitz::Array<double, 1>& atmosphere_bb);
  
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

  const boost::shared_ptr<Fo_Sswpgeometry_Master> geometry_interface() const { return geometry; }
  const boost::shared_ptr<Fo_Scalarss_Rtcalcs_Ilps> solar_interface() const { return solar_interface_; }

  bool do_deltam_scaling() const { return do_deltam_scaling_; }
  void do_deltam_scaling(bool do_scaling) { do_deltam_scaling_ = do_scaling; }
 
private:

  void init_interfaces(int nlayers, int surface_type);
  void copy_geometry_flags();

  const blitz::Array<double, 1> deltam_trunc_factor(const blitz::Array<double, 2>& pf) const;
  const blitz::Array<double, 2> deltam_linear_trunc_factor(const ArrayAd<double, 2>& pf) const;

  int num_layers_;
  int num_moments_;
  int num_streams_;
  int surface_type_;
  bool do_deltam_scaling_;

  blitz::Array<double, 1> height_diffs;

  boost::shared_ptr<Fo_Sswpgeometry_Master> geometry;
  boost::shared_ptr<Fo_Scalarss_Spherfuncs> legendre;
  boost::shared_ptr<Fo_Scalarss_Rtcalcs_Ilps> solar_interface_;

  // Specific implementation cast version of brdf_driver()
  boost::shared_ptr<LidortBrdfDriver> l_brdf_driver;
  FirstOrderDriver() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(FirstOrderDriver);
#endif
