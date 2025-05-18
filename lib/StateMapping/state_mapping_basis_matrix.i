%include "fp_common.i"

%{
#include "state_mapping_basis_matrix.h"
%}

%base_import(state_mapping)
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::StateMappingBasisMatrix);

namespace FullPhysics {
class StateMappingBasisMatrix : public StateMapping {
public:
  StateMappingBasisMatrix(const blitz::Array<double, 2> Basis_matrix);
  StateMappingBasisMatrix(const blitz::Array<double, 2> Basis_matrix,
			  const blitz::Array<double, 2>& Inverse_basis_matrix);
  virtual boost::shared_ptr<StateMapping> clone();
  %python_attribute(basis_matrix, blitz::Array<double, 2>);
  %python_attribute(inverse_basis_matrix, blitz::Array<double, 2>);
  %pickle_serialization();
};
}
