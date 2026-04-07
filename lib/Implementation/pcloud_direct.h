#ifndef PCLOUD_DIRECT_H
#define PCLOUD_DIRECT_H

#include "pcloud.h"
#include "sub_state_vector_array.h"
#include "double_with_unit.h"
#include "auto_derivative_with_unit.h"

namespace FullPhysics {
/****************************************************************//**
 Implements a direct representation of the Pcloud
 in the atmospheric state from the supplied value.
*******************************************************************/
class PcloudDirect :
    virtual public SubStateVectorArray<Pcloud> {
public:
  PcloudDirect(const ArrayWithUnit<double, 1>& pcloud, boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
  virtual ~PcloudDirect() {}

  virtual AutoDerivativeWithUnit<double>
  pressure_cloud(int sensor_index) const;
  virtual boost::shared_ptr<Pcloud> clone() const;

  virtual std::string sub_state_identifier() const
  { return "pcloud"; }

  std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const
  { Os << "PcloudDirect";}

private:
  Unit units;
  PcloudDirect() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
typedef SubStateVectorArray<Pcloud> SubStateVectorArrayPcloud;
}

FP_EXPORT_KEY(PcloudDirect)
FP_EXPORT_KEY(SubStateVectorArrayPcloud);
#endif
