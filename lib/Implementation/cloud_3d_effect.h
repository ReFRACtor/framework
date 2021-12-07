#ifndef CLOUD_3D_EFFECT_H
#define CLOUD_3D_EFFECT_H

#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include "spectrum_effect_imp_base.h"
#include "state_mapping.h"

namespace FullPhysics {

/****************************************************************//**
 Massie et al. (2020a) characterizes the biases resulting from 
 3D cloud effects resulting from clouds within about 10 km that 
 are not in the instrument’s field of view. The 3D-cloud effects 
 induce biases which are not corrected with current bias correction
 methods.

 This class uses the Schmidt et al. (2021) method to correct observed 
 radiances for these 3D effects by comparing 3D and 1D radiance 
 calculations. 
 
 The relationship between 3D and 1D radiances is:
 R_3d = R_1d(1 + cloud3d_offset + cloud3d_slope ∙ R_1d)

 where R is the reflectance, R = πI/F_"incident", I is the radiance, 
 and Fi incident is the incident irradiance, S*cos⁡(sza), where S
 is the solar spectrum.

 This class modifies the reflectance according to the above equation
 in a retrievable manner. The class is supplied with an initial guess
 for the offset and slope per band.

*******************************************************************/

class Cloud3dEffect : virtual public SpectrumEffectImpBase {
public:
  Cloud3dEffect(const double Offset, const double Slope, const std::string& Band_name, boost::shared_ptr<StateMapping> Mapping = boost::make_shared<StateMappingLinear>());

  virtual const AutoDerivative<double> offset() const;
  virtual const AutoDerivative<double> slope() const;

  virtual void apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& Forward_model_grid) const;

  virtual boost::shared_ptr<SpectrumEffect> clone() const;

  virtual std::string sub_state_identifier() const { return "cloud_3d/" + band_name; }

  virtual std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const;

  virtual std::string name() const { return "cloud_3d"; }

private:
  Cloud3dEffect() = default;
  std::string band_name;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(Cloud3dEffect);
#endif
