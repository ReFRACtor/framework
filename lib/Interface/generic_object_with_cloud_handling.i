// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "generic_object_with_cloud_handling.h"
%}

%base_import(generic_object)
%template(vector_generic_object_with_cloud_handling) std::vector<boost::shared_ptr<FullPhysics::GenericObjectWithCloudHandling> >;

%fp_shared_ptr(FullPhysics::GenericObjectWithCloudHandling);

namespace FullPhysics {

%feature("director") GenericObjectWithCloudHandling;
  
class GenericObjectWithCloudHandling : public GenericObject {
public:
  GenericObjectWithCloudHandling(bool Do_cloud = false);
  %python_attribute_with_set(do_cloud, bool);
  virtual void notify_do_cloud_update();
  virtual std::string desc() const;
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(generic_object_with_cloud_handling, GenericObjectWithCloudHandling)

// List of things "import *" will include
%python_export("GenericObjectWithCloudHandling");
