// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "ground_with_cloud_handling.h"
%}

%base_import(ground)
%base_import(generic_object_with_cloud_handling)

%fp_shared_ptr(FullPhysics::GroundWithCloudHandling);
namespace FullPhysics {
class GroundWithCloudHandling: public Ground, public Observer<Ground>,
  public GenericObjectWithCloudHandling {
public:
  GroundWithCloudHandling(const boost::shared_ptr<Ground> Ground_clear,
			  double Cloud_albedo, bool Do_cloud = false);
  virtual ArrayAd<double, 1> surface_parameter
    (const double wn, const int spec_index) const;
  %python_attribute_with_set(do_cloud, bool);
  %python_attribute(ground_clear, boost::shared_ptr<Ground>);
  %python_attribute_with_set(cloud_albedo, double);
};
}
