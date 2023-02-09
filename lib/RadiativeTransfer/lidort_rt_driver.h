#ifndef LIDORT_RT_DRIVER_H
#define LIDORT_RT_DRIVER_H

#include "multiscatt_rt_driver.h"
#include "lidort_brdf_driver.h"
#include "lidort_interface_masters.h"

namespace FullPhysics {

/****************************************************************//**
  LIDORT specific Radiative transfer interface implementation
 *******************************************************************/

class LidortRtDriver : public MultiScattRtDriver {
public:
  LidortRtDriver(int nstream, int nmoment, bool do_multi_scatt_only, int surface_type, 
          const blitz::Array<double, 1>& zen, bool pure_nadir, 
          bool do_solar_sources = true, bool do_thermal_emission = false, bool do_thermal_scattering = true);

  // Implements LIDORT specific extensions to MultiScattRtDriver actions
  void set_plane_parallel();
  void set_pseudo_spherical();
  void set_plane_parallel_plus_ss_correction();
  void set_line_of_sight();

  /// Abstract versions of interfaces
  const boost::shared_ptr<Spurr_Lps_Masters_Base> rt_interface() const 
  { return boost::dynamic_pointer_cast<Spurr_Lps_Masters_Base>(lidort_interface_); }

  const boost::shared_ptr<Spurr_Brdf_Lin_Sup_Masters_Base> brdf_interface() const
  { return boost::dynamic_pointer_cast<LidortBrdfDriver>(brdf_driver_)->brdf_interface(); }

  /// LIDORT specific interfaces
  const boost::shared_ptr<LidortBrdfDriver> lidort_brdf_driver() const
  { return boost::shared_ptr<LidortBrdfDriver>(boost::dynamic_pointer_cast<LidortBrdfDriver>(brdf_driver_)); }

  const boost::shared_ptr<Brdf_Lin_Sup_Masters> lidort_brdf_interface() const
  { return boost::dynamic_pointer_cast<LidortBrdfDriver>(brdf_driver_)->lidort_brdf_interface(); }

  const boost::shared_ptr<Lidort_Lps_Masters> lidort_interface() const { return lidort_interface_; }

  const blitz::Array<double, 1> get_intensity() const;
  void copy_jacobians(blitz::Array<double, 3>& jac_atm, blitz::Array<double, 2>& jac_surf_param, blitz::Array<double, 1>& jac_surf_temp, blitz::Array<double, 2>& jac_atm_temp) const;

protected:

  void initialize_rt(int nstream, int nmoment, bool do_solar_sources, bool do_thermal_emission, bool do_thermal_scattering);
  void copy_brdf_sup_outputs() const;
  void setup_phase_function(const blitz::Array<double, 3>& pf);
  void setup_linear_phase_function(const ArrayAd<double, 3>& pf);

  boost::shared_ptr<Lidort_Lps_Masters> lidort_interface_;

private:
  LidortRtDriver() {}
  void initialize_interface(int nstream, int nmoment);

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(LidortRtDriver);

#endif
