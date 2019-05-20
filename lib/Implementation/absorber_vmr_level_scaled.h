#ifndef ABSORBER_VMR_LEVEL_SCALED_H
#define ABSORBER_VMR_LEVEL_SCALED_H
#include "absorber_vmr_level.h"

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This particular implementation uses a passed vmr profile
  (interpolated to the current pressure grid), along with a scale factor.
*******************************************************************/
class AbsorberVmrLevelScaled : public AbsorberVmrLevel {
public:
  AbsorberVmrLevelScaled(const boost::shared_ptr<Pressure>& Press,
                         const blitz::Array<double, 1>& Vmr_profile,
                         double Scale,                         
                         bool Scale_flag,
                         const std::string& Gas_name);
  virtual ~AbsorberVmrLevelScaled() {}
};
}
#endif
