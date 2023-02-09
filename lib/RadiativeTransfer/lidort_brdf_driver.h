#ifndef LIDORT_BRDF_DRIVER_H
#define LIDORT_BRDF_DRIVER_H

#include "multiscatt_brdf_driver.h"
#include "lidort_interface_masters.h"

namespace FullPhysics {

/****************************************************************//**
  LIDORT specific BRDF driver implementation
 *******************************************************************/

class LidortBrdfDriver : public MultiScattBrdfDriver {
public:
  LidortBrdfDriver(int nstream, int nmoment);
  virtual ~LidortBrdfDriver() = default;

  /// Interface to the underlying BRDF interface, implements required interface
  const boost::shared_ptr<Spurr_Brdf_Lin_Sup_Masters_Base> brdf_interface() const { return brdf_interface_; }

  /// Interface to the underlying LIDORT specific BRDF interface
  const boost::shared_ptr<Brdf_Lin_Sup_Masters> lidort_brdf_interface() const { return brdf_interface_; }

private:
  LidortBrdfDriver() = default;

  boost::shared_ptr<Brdf_Lin_Sup_Masters> brdf_interface_;
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(LidortBrdfDriver);

#endif
