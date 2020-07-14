#ifndef OPTICAL_PROP_INIT_BASE_H
#define OPTICAL_PROP_INIT_BASE_H

#include "optical_properties_imp_base.h"

#include "absorber.h"
#include "rayleigh.h"
#include "aerosol.h"

namespace FullPhysics {

/****************************************************************//**
 Extension of OpticalPropertiesImpBase which allows inheriting
 classes to share common value extraction initialization methods.
 *******************************************************************/

class OpticalPropertiesInitBase : public virtual OpticalPropertiesImpBase {
public:
  OpticalPropertiesInitBase() : OpticalPropertiesImpBase() { }

  virtual ~OpticalPropertiesInitBase() = default;

  virtual void initialize(const ArrayAd<double, 1>& rayleigh_od, 
			  const ArrayAd<double, 2>& gas_od);

  virtual void initialize(const ArrayAd<double, 1>& rayleigh_od, 
			  const ArrayAd<double, 2>& gas_od,
			  const ArrayAd<double, 2>& aerosol_ext_od,
			  const ArrayAd<double, 2>& aerosol_sca_od,
			  const std::vector<ArrayAd<double, 3> >& aerosol_pf_moments);

  virtual void initialize(const DoubleWithUnit spectral_point,
			  const int channel_index,
			  const boost::shared_ptr<Absorber>& absorber,
			  const boost::shared_ptr<Rayleigh>& rayleigh,
			  const boost::shared_ptr<Aerosol>& aerosol);

protected:

  virtual void initialize_with_jacobians(const ArrayAd<double, 1>& rayleigh_od, 
					 const ArrayAd<double, 2>& gas_od,
					 const ArrayAd<double, 2>& aerosol_ext_od,
					 const ArrayAd<double, 2>& aerosol_sca_od,
					 const boost::shared_ptr<AerosolPhaseFunctionHelper>& aer_pf_helper) = 0;
private:  
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(OpticalPropertiesInitBase);
#endif
