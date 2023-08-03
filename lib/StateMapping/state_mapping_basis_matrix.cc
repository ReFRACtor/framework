#include "state_mapping_basis_matrix.h"
#include "linear_algebra.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StateMappingBasisMatrix::serialize(Archive& ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateMapping)
    & FP_NVP_(basis_matrix) & FP_NVP_(inverse_basis_matrix);
}

FP_IMPLEMENT(StateMappingBasisMatrix);
#endif

//-----------------------------------------------------------------------
/// Constructor. This takes a matrix that maps from the retrieval grid
/// to the forward model grid. This will be m x n, with m the size of 
/// the forward grid and n the size of the retrieval grid (so m > n).
///
/// When comparing with py_retrieve, note that this is the transpose
/// of the basis matrix it uses.
//-----------------------------------------------------------------------

StateMappingBasisMatrix::StateMappingBasisMatrix
(const blitz::Array<double, 2> Basis_matrix)
  :  basis_matrix_(Basis_matrix.copy())
{
  if(Basis_matrix.rows() == 0 || Basis_matrix.cols() ==0)
    throw Exception("Basis matrix is size 0");
  inverse_basis_matrix_.reference(generalized_inverse(basis_matrix_));
}

//-----------------------------------------------------------------------
/// Go to forward model grid
//-----------------------------------------------------------------------

ArrayAd<double, 1>
StateMappingBasisMatrix::mapped_state(const ArrayAd<double, 1>& retrieval_values) const
{
  blitz::firstIndex i1;
  blitz::secondIndex i2;
  blitz::Array<AutoDerivative<double>, 1> r(basis_matrix_.rows());
  blitz::Array<AutoDerivative<double>, 2> a(basis_matrix_.shape());
  a = auto_derivative(basis_matrix_);
  blitz::Array<AutoDerivative<double>, 1> b = retrieval_values.to_array();
  r = blitz::sum(a(i1,i2) * b(i2), i2);
  return r;
}

//-----------------------------------------------------------------------
/// Go to retrieval grid
//-----------------------------------------------------------------------

ArrayAd<double, 1> StateMappingBasisMatrix::retrieval_state
(const ArrayAd<double, 1>& initial_values) const
{
  blitz::firstIndex i1;
  blitz::secondIndex i2;
  blitz::Array<AutoDerivative<double>, 1> r(inverse_basis_matrix_.rows());
  blitz::Array<AutoDerivative<double>, 2> a(inverse_basis_matrix_.shape());
  a = auto_derivative(inverse_basis_matrix_);
  blitz::Array<AutoDerivative<double>, 1> b = initial_values.to_array();
  r = blitz::sum(a(i1,i2) * b(i2), i2);
  return r;
}
