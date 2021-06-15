// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "generic_object_with_cloud_handling.h"
%}

%base_import(generic_object)
%fp_shared_ptr(FullPhysics::GenericObjectWithCloudHandling);

namespace FullPhysics {

%feature("director") GenericObjectWithCloudHandling;
  
class GenericObjectWithCloudHandling : public GenericObject {
public:
  GenericObjectWithCloudHandling(bool Do_cloud);
  %python_attribute_with_set(do_cloud, bool);
  virtual void notify_do_cloud_update();
  std::string print_to_string() const;
  %pickle_serialization();
};
}

%template(vector_generic_object_with_cloud_handling) std::vector<boost::shared_ptr<FullPhysics::GenericObjectWithCloudHandling> >;
