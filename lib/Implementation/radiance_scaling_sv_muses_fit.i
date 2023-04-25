// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "radiance_scaling_sv_muses_fit.h"
#include "sub_state_vector_array.h"
%}
%base_import(radiance_scaling)
%base_import(instrument_correction)
%base_import(sub_state_vector_array)
%import "double_with_unit.i"
%import "auto_derivative.i"
%import "array_ad.i"
%import "spectral_domain.i"
%import "spectral_range.i"

%fp_shared_ptr(FullPhysics::RadianceScalingSvMusesFit);

namespace FullPhysics {

// Force to be not abstract, SWIG had troubles seeing that the clone methods ARE implemented below
%feature("notabstract") RadianceScalingSvMusesFit;

class RadianceScalingSvMusesFit : public RadianceScaling, 
                             public SubStateVectorArray<InstrumentCorrection> {
public:
  RadianceScalingSvMusesFit(const blitz::Array<double, 1>& Coeff, 
                       const DoubleWithUnit& Band_ref,
                       const std::string& Band_name);
  %python_attribute(sub_state_identifier, std::string);
  virtual std::string state_vector_name_i(int i) const;
  virtual boost::shared_ptr<InstrumentCorrection> clone() const;
  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const;
  virtual void print(std::ostream& Os) const;
  virtual void notify_update(const StateVector& Sv);
  %python_attribute(radiance_scaling_coeff_uncertainty, blitz::Array<double, 1>)
  %pickle_serialization();
};
}

