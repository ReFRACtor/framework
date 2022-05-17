// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "pressure_with_cloud_handling.h"
%}

%base_import(pressure)
%base_import(generic_object_with_cloud_handling)
%import "auto_derivative.i"
%fp_shared_ptr(FullPhysics::PressureWithCloudHandling);

namespace FullPhysics {
class PressureWithCloudHandling : public Pressure, public Observer<Pressure>,
  public GenericObjectWithCloudHandling {
public:
  PressureWithCloudHandling(const boost::shared_ptr<Pressure> Press_clear,
                            double Cloud_pressure_level, bool do_cloud = false);
  %python_attribute(pressure_clear, boost::shared_ptr<Pressure>);
  %python_attribute_derived(type_preference, TypePreference)
  virtual ArrayAdWithUnit<double, 1>
  pressure_grid(PressureGridType Gtype = INCREASING_PRESSURE) const;
  virtual boost::shared_ptr<Pressure> clone() const;
  %python_attribute_with_set(cloud_pressure_level, double);
  %pickle_serialization();
};
}
