#ifndef SPURR_RT_DRIVER_H
#define SPURR_RT_DRIVER_H

#include <boost/shared_ptr.hpp>
#include <blitz/array.h>
#include "array_ad.h"
#include "generic_object.h"
#include "rt_atmosphere.h"

#include "spurr_brdf_driver.h"

/****************************************************************//**
  Contains classes to abstract away details in various Spurr
  Radiative Transfer software.
*******************************************************************/

namespace FullPhysics {

/****************************************************************//**
  Abstracts away set up of Radiative Transfer software from Rob
  Spurr into a simpler common inteface used by the L2 software

  This interface should be independent of the L2 Atmosphere
  class to make unit testing easier.
*******************************************************************/

class SpurrRtDriver : public Printable<SpurrRtDriver> {

public:

  SpurrRtDriver(bool do_solar = true, bool do_thermal = false) 
      : do_solar_sources(do_solar), do_thermal_emission(do_thermal) {}

  //-----------------------------------------------------------------------
  /// Some drivers may have assumptions about the number of layers or
  /// surface type. They can override this function to handle any
  /// changes to the atm. Default is to do nothing.
  //-----------------------------------------------------------------------
  
  virtual void notify_update(const RtAtmosphere& UNUSED(atm)) { }

  /// Computes reflectance without jacobians
  virtual const blitz::Array<double, 1> reflectance_calculate(const blitz::Array<double, 1>& height_grid,
                                                              double sza, double azm, double zen,
                                                              int surface_type,
                                                              const blitz::Array<double, 1>& surface_parameters,
                                                              const blitz::Array<double, 1>& od, 
                                                              const blitz::Array<double, 1>& ssa,
                                                              const blitz::Array<double, 3>& pf,
                                                              double surface_bb = 0,
                                                              const blitz::Array<double, 1>& atmosphere_bb = blitz::Array<double,1>());
  
  // Computes reflectance and jacobians for profiles as well as surface
  virtual void reflectance_and_jacobian_calculate(const blitz::Array<double, 1>& height_grid,
                                                  double sza, double azm, double zen,
                                                  int surface_type,
                                                  ArrayAd<double, 1>& surface_parameters,
                                                  const ArrayAd<double, 1>& od, 
                                                  const ArrayAd<double, 1>& ssa,
                                                  const ArrayAd<double, 3>& pf,
                                                  blitz::Array<double, 1>& reflectance,
                                                  blitz::Array<double, 3>& jac_atm, 
                                                  blitz::Array<double, 2>& jac_surf_param,
                                                  blitz::Array<double, 1>& jac_surf_temp,
                                                  blitz::Array<double, 2>& jac_atm_temp,
                                                  double surface_bb = 0,
                                                  const blitz::Array<double, 1>& atmosphere_bb = blitz::Array<double,1>());

  /// Access to BRDF driver
  const boost::shared_ptr<SpurrBrdfDriver> brdf_driver() const { return brdf_driver_; }

  /// Setup height grid, should only be called once per instance or if
  /// the height grid changes
  virtual void setup_height_grid(const blitz::Array<double, 1>& height_grid)  = 0;

  /// Setup viewing geometry, should only be called once per instance or if
  /// the viewing geometry changes
  virtual void setup_geometry(double sza, double azm, double zen) = 0;

  /// Set up thermal emission inputs
  virtual void setup_thermal_inputs(double surface_bb, const blitz::Array<double, 1>& atmosphere_bb) = 0;

  /// Set up optical depth, single scattering albedo and phase function
  /// Should be called per spectral point
  virtual void setup_optical_inputs(const blitz::Array<double, 1>& od, 
                                    const blitz::Array<double, 1>& ssa,
                                    const blitz::Array<double, 3>& pf) = 0;
  
  /// Mark that we are not retrieving weighting functions
  virtual void clear_linear_inputs() =  0;

  /// Set up linearization, weighting functions
  virtual void setup_linear_inputs(const ArrayAd<double, 1>& od,
                                   const ArrayAd<double, 1>& ssa,
                                   const ArrayAd<double, 3>& pf,
                                   bool do_surface_linearization) = 0;

  /// Perform radiative transfer calculation with the values
  /// setup by setup_optical_inputs and setup_linear_inputs
  virtual void calculate_rt() const = 0;

  /// Retrieve the intensity value calculated
  virtual const blitz::Array<double, 1> get_intensity() const = 0;

  /// Copy jacobians out of internal xdata structures
  virtual void copy_jacobians(blitz::Array<double, 3>& jac_atm, blitz::Array<double, 2>& jac_surf_params, blitz::Array<double, 1>& jac_surf_temp, blitz::Array<double, 2>& jac_atm_temp) const = 0;
  virtual void print(std::ostream& Os) const {Os << "SpurrRtDriver";}

protected:

  bool do_solar_sources, do_thermal_emission;

  /// Spurr BRDF class interface class to use
  mutable boost::shared_ptr<SpurrBrdfDriver> brdf_driver_;
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(SpurrRtDriver);

#endif
