%include "fp_common.i"

%{
#include "spectral_window_range.h"
%}

%base_import(spectral_window)

%import "array_with_unit.i"
%import "double_with_unit.i"
%import "sample_grid.i"

%fp_shared_ptr(FullPhysics::SpectralWindowRange);

namespace FullPhysics {

class SampleGrid;

class SpectralWindowRange : public SpectralWindow {
public:
  SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges);
  SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges, const blitz::Array<bool, 2>& Bad_sample_mask);
  SpectralWindowRange(const ArrayWithUnit<double, 3>& Microwindow_ranges, const std::vector<blitz::Array<bool, 1> >& Bad_sample_mask);
  virtual std::vector<int> grid_indexes(const SpectralDomain& Grid, int Spec_index) const;
  %python_attribute_derived(number_spectrometer, int);
  %python_attribute_derived(spectral_bound, SpectralBound);
  %python_attribute_with_set(range_array, ArrayWithUnit<double, 3>);
  const blitz::Array<bool, 1>& bad_sample_mask(int sensor_index) const;
  void bad_sample_mask(const blitz::Array<bool, 1>& mask, int sensor_index);
  %python_attribute_with_set(dispersion, std::vector<boost::shared_ptr<SampleGrid> >)
  %pickle_serialization();
};

}

// List of things "import *" will include
%python_export("SpectralWindowRange");
