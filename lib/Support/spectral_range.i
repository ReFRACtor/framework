// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "spectral_range.h"
%}
%base_import(generic_object)
%import "array_with_unit.i"
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::SpectralRange)
namespace FullPhysics {
class SpectralRange : public GenericObject {
public:
  SpectralRange(const ArrayWithUnit<double, 1>& Data);
  SpectralRange(const blitz::Array<double, 1>& FORCE_COPY, 
		const Unit& U);
  SpectralRange(const blitz::Array<double, 1>& FORCE_COPY, 
		const std::string& U);
  SpectralRange(const ArrayAd<double, 1>& Data, 
		const Unit& U);
  SpectralRange(const ArrayAd<double, 1>& Data, 
		const std::string& U);
  SpectralRange(const blitz::Array<double, 1>& FORCE_COPY, 
		const Unit& U,
		const blitz::Array<double, 1>& FORCE_COPY);
  SpectralRange(const blitz::Array<double, 1>& FORCE_COPY, 
		const std::string& U,
		const blitz::Array<double, 1>& FORCE_COPY);
  SpectralRange(const ArrayAd<double, 1>& Data, 
		const Unit& U,
		const blitz::Array<double, 1>& FORCE_COPY);
  %python_attribute(data, blitz::Array<double, 1>)
  %python_attribute(uncertainty, blitz::Array<double, 1>)
  %python_attribute(units, Unit)
  %python_attribute(data_ad, ArrayAd<double, 1>)
  SpectralRange convert(const Unit& R) const;
  SpectralRange convert(const std::string& R) const;
  std::string print_to_string() const;
  %pickle_serialization();

  %pythoncode {
def copy(self):
    return self.__class__(self.data_ad.copy(), self.units)
  }

};
}

%template(vector_SpectralRange) std::vector<boost::shared_ptr<FullPhysics::SpectralRange> >;
