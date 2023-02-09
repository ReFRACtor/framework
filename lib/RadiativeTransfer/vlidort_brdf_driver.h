#ifndef VLIDORT_BRDF_DRIVER_H
#define VLIDORT_BRDF_DRIVER_H

#include "multiscatt_brdf_driver.h"
#include "vlidort_interface_masters.h"

namespace FullPhysics {

/****************************************************************//**
  VLIDORT specific BRDF driver implementation
 *******************************************************************/

class VLidortBrdfDriver : public MultiScattBrdfDriver {
public:
  VLidortBrdfDriver(int nstream, int nmoment);
  virtual ~VLidortBrdfDriver() = default;

  /// Interface to the underlying BRDF interface, implements required interface
  const boost::shared_ptr<Spurr_Brdf_Lin_Sup_Masters_Base> brdf_interface() const { return brdf_interface_; }

  /// Interface to the underlying VLIDORT specific BRDF interface
  const boost::shared_ptr<VBrdf_Linsup_Masters> vlidort_brdf_interface() const { return brdf_interface_; };

private:
  VLidortBrdfDriver() = default;

  boost::shared_ptr<VBrdf_Linsup_Masters> brdf_interface_;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

}

FP_EXPORT_KEY(VLidortBrdfDriver);

#endif
