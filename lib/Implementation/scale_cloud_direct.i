
%include "fp_common.i"

%{
#include "scale_cloud_direct.h"
%}

%base_import(scale_cloud)
%import "sub_state_vector_array.i"

%fp_shared_ptr(FullPhysics::ScaleCloudDirect);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::ScaleCloud>)

%template(SubStateVectorArrayScaleCloudDirect) FullPhysics::SubStateVectorArray<FullPhysics::ScaleCloud>;

namespace FullPhysics {
class ScaleCloudDirect : public SubStateVectorArray<ScaleCloud> {
public:
  ScaleCloudDirect(const blitz::Array<double, 1>& scale_cloud);
  ScaleCloudDirect(const blitz::Array<double, 1>& scale_cloud, boost::shared_ptr<StateMapping> in_map);
  virtual ~ScaleCloudDirect() {}
  virtual AutoDerivative<double>
    scale_cloud(int sensor_index) const;
  virtual boost::shared_ptr<ScaleCloud> clone() const;

  std::string state_vector_name_i(int i) const;
  %pickle_serialization();
};
}
