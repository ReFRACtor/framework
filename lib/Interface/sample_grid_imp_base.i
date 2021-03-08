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

  // In Python you need to override the renamed function _v_sample_grid
  %python_attribute_abstract(sample_grid, SpectralDomain);

  %sub_state_virtual_func(SampleGrid);
  %pickle_serialization();
protected:
  SampleGridImpBase(const blitz::Array<double, 1>& Coeff,
                    boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
};

}
