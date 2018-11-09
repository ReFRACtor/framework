%include <std_vector.i>
%include "common.i"

%{
#include "ils.h"
%}

%base_import(state_vector_observer)

%import "array_ad.i"
%import "spectral_domain.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::Ils);

namespace FullPhysics {

%feature("director") Ils;

class Ils : public StateVectorObserver {
public:
  virtual ~Ils();
  std::string print_to_string() const;
  virtual blitz::Array<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const blitz::Array<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const = 0;
  virtual ArrayAd<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const ArrayAd<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const = 0;
  virtual boost::shared_ptr<Ils> clone() const = 0;

  // Needed for directors, can not use %python_attribute here or else there will
  // be missing symbol problems in the director
  virtual const SpectralDomain pixel_grid() const = 0;
  virtual const DoubleWithUnit ils_half_width() const = 0;
  virtual void ils_half_width(const DoubleWithUnit& half_width) = 0;
};
}

%template(vector_ils) std::vector<boost::shared_ptr<FullPhysics::Ils> >;
