
%include "fp_common.i"

%{
#include "pcloud_direct.h"
%}

%base_import(pcloud)
%import "sub_state_vector_array.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::PcloudDirect);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::Pcloud>)

%template(SubStateVectorArrayPcloudDirect) FullPhysics::SubStateVectorArray<FullPhysics::Pcloud>;

namespace FullPhysics {
class PcloudDirect : public SubStateVectorArray<Pcloud> {
public:
  PcloudDirect(const ArrayWithUnit<double, 1>& pcloud);
  PcloudDirect(const ArrayWithUnit<double, 1>& pcloud, boost::shared_ptr<StateMapping> in_map);
  virtual ~PcloudDirect() {}
  virtual AutoDerivativeWithUnit<double>
    pressure_cloud(int sensor_index) const;
  virtual boost::shared_ptr<Pcloud> clone() const;

  std::string state_vector_name_i(int i) const;
  %pickle_serialization();
};
}

// List of things "import *" will include
%python_export("PcloudDirect");
