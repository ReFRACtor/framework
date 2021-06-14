// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "pressure_with_cloud_handling.h"
%}

%base_import(pressure)
%import "auto_derivative.i"
%fp_shared_ptr(FullPhysics::PressureWithCloudHandling);

namespace FullPhysics {
class PressureWithCloudHandling : public Pressure, public Observer<Pressure> {
public:
  PressureWithCloudHandling(const boost::shared_ptr<Pressure> Press_clear,
			    double Cloud_pressure_level, bool do_cloud = false);
  %python_attribute_with_set(do_cloud, bool);
  %python_attribute(pressure_clear, boost::shared_ptr<Pressure>);
  %python_attribute_derived(pressure_grid, ArrayAdWithUnit<double, 1>);
  virtual boost::shared_ptr<Pressure> clone() const;
  %python_attribute_with_set(cloud_pressure_level, double);
  %pickle_serialization();
};
}
