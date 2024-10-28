%include "fp_common.i"

%{
#include "sample_grid_imp_base.h"
%}

%base_import(state_vector)
%base_import(sub_state_vector_array)
%base_import(sample_grid)
%base_import(state_mapping)
%base_import(state_mapping_linear)

%fp_shared_ptr(FullPhysics::SampleGridImpBase);
%fp_shared_ptr(FullPhysics::SubStateVectorArray<FullPhysics::SampleGrid>);
%nodefaultctor FullPhysics::SubStateVectorArray<FullPhysics::SampleGrid>;

namespace FullPhysics {

%template(SubStateVectorArraySampleGrid) FullPhysics::SubStateVectorArray<SampleGrid>;

// Allow these classes to be derived from in Python.
%feature("director") SampleGridImpBase;

class SampleGridImpBase: public SubStateVectorArray<SampleGrid> {
public:
  virtual boost::shared_ptr<SampleGrid> clone() const = 0;
  %sub_state_virtual_func(SampleGrid);
  // In Python you need to override the renamed function _v_sample_grid
  %python_attribute_abstract(sample_grid, SpectralDomain);
  %python_attribute(pixel_grid, SpectralDomain);
  %python_attribute_with_set(sv_name, std::vector<std::string>);
  virtual std::string desc() const;
  std::string print_to_string() const;
  std::string print_parent() const;
  %pickle_serialization();
protected:
  SampleGridImpBase();
  SampleGridImpBase(const blitz::Array<double, 1>& Coeff,
                    boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
};

}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(sample_grid_imp_base, SampleGridImpBase)

// List of things "import *" will include
%python_export("SampleGridImpBase",  "SubStateVectorArraySampleGrid");

