#ifndef ABSORBER_VMR_LEVEL_LOG_H
#define ABSORBER_VMR_LEVEL_LOG_H

#include "absorber_vmr_level.h"
#include <boost/lexical_cast.hpp>

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This particular implementation represents VMR values in log space
  in the state vector.
*******************************************************************/
class AbsorberVmrLevelLog : virtual public AbsorberVmrLevel {
public:
  AbsorberVmrLevelLog(const boost::shared_ptr<Pressure>& Press,
		      const blitz::Array<double, 1>& Vmr,
		      const blitz::Array<bool, 1>& Vmr_flag,
		      const std::string& Gas_name);

  virtual ~AbsorberVmrLevelLog() = default;
    
  virtual boost::shared_ptr<AbsorberVmr> clone() const;
protected:
  AbsorberVmrLevelLog() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(AbsorberVmrLevelLog);
#endif
