#ifndef MULTISCATT_RT_DRIVER_H
#define MULTISCATT_RT_DRIVER_H

#include "spurr_rt_driver.h"
#include "spurr_interface_types.h"
#include "spurr_interface_masters.h"

/****************************************************************//**
  Contains classes to abstract away details of the common behaviors
  for setting up LIDORT and VLIDORT multiple scattering code
*******************************************************************/

namespace FullPhysics {

/****************************************************************//**
  Abstracts away set up of LIDORT and VLIDOR Radiative Transfer 
  software from Rob Spurr into a simpler common inteface to 
  reduce duplication of the setup to both multiple scattering codes.
*******************************************************************/

class MultiScattRtDriver : public SpurrRtDriver {

public:

  MultiScattRtDriver(bool do_solar = true, bool do_thermal = false) :
    SpurrRtDriver(do_solar, do_thermal) {}

  virtual int number_moment() const;
  virtual int number_stream() const;
  virtual int surface_type() const { return surface_type_; }
  virtual bool do_multi_scatt_only() const { return do_multi_scatt_only_; }
  virtual bool pure_nadir() const { return pure_nadir_; }

  void setup_sphericity(blitz::Array<double, 1> zen, bool do_multi_scatt_only, bool pure_nadir);
  void set_plane_parallel();
  void set_pseudo_spherical();
  void set_plane_parallel_plus_ss_correction();
  void set_line_of_sight();

  virtual void notify_update(const RtAtmosphere& atm);

  virtual const boost::shared_ptr<Spurr_Lps_Masters_Base> rt_interface() const = 0;
  virtual const boost::shared_ptr<Spurr_Brdf_Lin_Sup_Masters_Base> brdf_interface() const = 0;

protected:

  virtual void initialize_interface(int nstream, int nmoment) = 0;

  virtual void initialize_brdf(int surface_type);
  virtual void initialize_rt(int nstream, int nmoment, bool do_solar_sources, bool do_thermal_emission, bool do_thermal_scattering);

  virtual void copy_brdf_sup_outputs() const = 0;

private:

  int surface_type_;
  bool do_multi_scatt_only_;
  bool pure_nadir_;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(MultiScattRtDriver);

#endif
