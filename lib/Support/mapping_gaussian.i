%include "fp_common.i"

%{
#include "mapping_gaussian.h"
%}

%base_import(mapping)
%import "array_ad.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::MappingGaussian);


namespace FullPhysics {

class MappingGaussian : public Mapping  {
public:
  MappingGaussian(const boost::shared_ptr<Pressure>& in_press, bool Linear_AOD,
		  double Min_Desired = 1e-9);
  virtual boost::shared_ptr<Mapping> clone();
  %python_attribute(is_linear_total, bool);
  %python_attribute(min_desired, double);
  %python_attribute(pressure, boost::shared_ptr<Pressure>);
  virtual AutoDerivative<double>
    total_optical_depth(const ArrayAd<double, 1>& component) const;
  %pickle_serialization();
};
}

