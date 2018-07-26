// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "sample_grid.h"
%}
%base_import(generic_object)
%import "spectral_domain.i"

%fp_shared_ptr(FullPhysics::SampleGrid);

namespace FullPhysics {
class SampleGrid : public GenericObject {
public:
  virtual ~SampleGrid();
  virtual SpectralDomain sample_grid(int spec_index) const = 0;
};
}
