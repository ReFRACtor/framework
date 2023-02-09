#ifndef MULTISCATT_BRDF_DRIVER_H
#define MULTISCATT_BRDF_DRIVER_H

#include "spurr_brdf_driver.h"
#include "spurr_interface_masters.h"

namespace FullPhysics {

/****************************************************************//**
  Common behavior of LIDORT and VLIDORT BRDF drivers
 *******************************************************************/

class MultiScattBrdfDriver : public SpurrBrdfDriver {
public:

  virtual ~MultiScattBrdfDriver() = default;

  /// Interface to the underlying BRDF interface
  virtual const boost::shared_ptr<Spurr_Brdf_Lin_Sup_Masters_Base> brdf_interface() const = 0;

  virtual void setup_geometry(double sza, double azm, double zen);

  virtual int n_brdf_kernels() const;

  virtual int n_kernel_factor_wfs() const;
  virtual int n_kernel_params_wfs() const;
  virtual int n_surface_wfs() const;
  virtual bool do_kparams_derivs(const int kernel_index) const;
  virtual bool do_shadow_effect() const;

protected:
  MultiScattBrdfDriver(int nstream, int nmoment);
  MultiScattBrdfDriver() = default;

  int nstream_;
  int nmoment_;

  void init();

  virtual void calculate_brdf() const;

  virtual void n_brdf_kernels(const int n_kernels);
  virtual void n_kernel_factor_wfs(const int n_factors);
  virtual void n_kernel_params_wfs(const int n_params);
  virtual void n_surface_wfs(const int n_wfs);
  virtual void do_kparams_derivs(const int kernel_index, const bool do_kparams);
  virtual void do_shadow_effect(const bool do_shadow) const;

  virtual void initialize_kernel_parameters(const int kernel_index,
                                            const int which_brdf,
                                            const bool lambertian_flag,
                                            const int n_brdf_parameters,
                                            const bool do_factor_wfs,
                                            const blitz::Array<bool, 1>& do_params_wfs);

private:

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(MultiScattBrdfDriver);

#endif
