// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "empirical_orthogonal_function.h"
%}
%base_import(sub_state_vector_array)
%base_import(instrument_correction)
%import "array_with_unit.i"
%import "sample_grid.i"
%import "hdf_file.i"
%fp_shared_ptr(FullPhysics::EmpiricalOrthogonalFunction);

namespace FullPhysics {
class EmpiricalOrthogonalFunction: 
  public SubStateVectorArray<InstrumentCorrection> {
public:
  EmpiricalOrthogonalFunction(double Coeff, 
                              const ArrayWithUnit<double, 1>& Eof_waveform,
                              int Order,
                              const std::string& Band_name,
                              const std::string& Hdf_group = "N/A");
  EmpiricalOrthogonalFunction(double Coeff, 
                              const SampleGrid& Disp,
                              const HdfFile& Hdf_static_input,
                              int Spec_index,
                              int Sounding_number,
                              int Order,
                              const std::string& Band_name,
                              const std::string& Hdf_group = 
                              "Instrument/EmpiricalOrthogonalFunction_1");
  EmpiricalOrthogonalFunction(double Coeff, 
                              const HdfFile& Hdf_static_input,
                              int Spec_index,
                              int Sounding_number,
                              int Order,
                              const std::string& Band_name,
                              const std::string& Hdf_group = 
                              "Instrument/EmpiricalOrthogonalFunction_1");
  EmpiricalOrthogonalFunction(double Coeff, 
                              const HdfFile& Hdf_static_input,
                              const ArrayWithUnit<double, 1>& Uncertainty,
                              int Spec_index,
                              int Sounding_number,
                              int Order,
                              const std::string& Band_name,
                              const std::string& Hdf_group = 
                              "Instrument/EmpiricalOrthogonalFunction",
                              double Scale_to_stddev = 1e19);
  virtual boost::shared_ptr<InstrumentCorrection> clone() const;
  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const;
  %python_attribute(eof, ArrayWithUnit<double, 1>)
  %python_attribute(order, int)
  %sub_state_virtual_func(InstrumentCorrection);
};
}
