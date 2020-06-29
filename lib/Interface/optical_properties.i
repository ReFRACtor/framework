// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "optical_properties.h"
#include "sub_state_vector_array.h"
%}

%base_import(generic_object)
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::OpticalProperties)

namespace FullPhysics {
class OpticalProperties : public GenericObject {
public:
  virtual int number_layers() const = 0;
  virtual int number_gas_particles() const = 0;
  virtual int number_aerosol_particles() const = 0;
  virtual ArrayAd<double, 1> rayleigh_optical_depth() const = 0;
  virtual ArrayAd<double, 2> gas_optical_depth_per_particle() const = 0;
  virtual ArrayAd<double, 2> aerosol_extinction_optical_depth_per_particle() const = 0;
  virtual ArrayAd<double, 2> aerosol_scattering_optical_depth_per_particle() const = 0;
  virtual const std::vector<ArrayAd<double, 3> > aerosol_phase_function_moments_per_particle(int num_moments = -1, int num_scattering = -1) const = 0;
  virtual ArrayAd<double, 1> gas_optical_depth_per_layer() const = 0;
  virtual ArrayAd<double, 1> aerosol_extinction_optical_depth_per_layer() const = 0;
  virtual ArrayAd<double, 1> aerosol_scattering_optical_depth_per_layer() const = 0;
  virtual ArrayAd<double, 1> total_optical_depth() const = 0;
  virtual ArrayAd<double, 1> total_single_scattering_albedo() const = 0;
  virtual ArrayAd<double, 1> rayleigh_fraction() const = 0;
  virtual ArrayAd<double, 2> aerosol_fraction() const = 0;
  virtual ArrayAd<double, 3> rayleigh_phase_function_moments_portion(int num_scattering) const = 0;
  virtual ArrayAd<double, 3> aerosol_phase_function_moments_portion(int num_moments = -1, int num_scattering = -1) const = 0;
  virtual ArrayAd<double, 3> total_phase_function_moments(int num_moments = -1, int num_scattering = -1) const = 0;
  virtual blitz::Array<double, 3> intermediate_jacobian() const = 0;
  std::string print_to_string() const;
  %pickle_serialization();
};
}
  
