#ifndef ABSORBER_VMR_LEVEL_SCALED2_H
#define ABSORBER_VMR_LEVEL_SCALED2_H

#include "absorber_vmr_level.h"
#include <boost/lexical_cast.hpp>

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This particular implementation uses a passed vmr profile
  (interpolated to the current pressure grid), along with a scale factor.
*******************************************************************/
class AbsorberVmrLevelScaled2 : public AbsorberVmrLevel {
public:
  AbsorberVmrLevelScaled2(const boost::shared_ptr<Pressure>& Press,
                         const blitz::Array<double, 1>& Vmr_profile,
                         double Scale,
                         const blitz::Array<bool, 1>& Scale_flag,
                         const std::string& Gas_name);

    virtual ~AbsorberVmrLevelScaled2() = default;
};
}
#endif
