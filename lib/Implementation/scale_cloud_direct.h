#ifndef SCALE_CLOUD_DIRECT_H
#define SCALE_CLOUD_DIRECT_H

#include "scale_cloud.h"
#include "sub_state_vector_array.h"
#include "auto_derivative.h"

namespace FullPhysics {
/****************************************************************//**
 Implements a direct representation of the surface temperature
 in the atmospheric state from the supplied value.
*******************************************************************/
class ScaleCloudDirect :
    virtual public SubStateVectorArray<ScaleCloud> {
public:
  ScaleCloudDirect(const blitz::Array<double, 1>& scale_cloud, boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
  virtual ~ScaleCloudDirect() {}

  virtual AutoDerivative<double>
  scale_cloud(int sensor_index) const;
  virtual boost::shared_ptr<ScaleCloud> clone() const;

  virtual std::string sub_state_identifier() const
  { return "scale_cloud"; }

  std::string state_vector_name_i(int i) const;

  virtual void print(std::ostream& Os) const
  { Os << "ScaleCloudDirect";}
private:
  ScaleCloudDirect() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
typedef SubStateVectorArray<ScaleCloud> SubStateVectorArrayScaleCloud;
}

FP_EXPORT_KEY(ScaleCloudDirect)
FP_EXPORT_KEY(SubStateVectorArrayScaleCloud);
#endif
