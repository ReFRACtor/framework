#ifndef STATE_MAPPING_BASIS_MATRIX_H
#define STATE_MAPPING_BASIS_MATRIX_H

#include <blitz/array.h>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "array_ad.h"
#include "state_mapping.h"


namespace FullPhysics {
/****************************************************************//**
  This class implements a mapping going from a retrieval to a forward
  model grid using a basis matrix. This is what muses_py does - this
  is somewhat like and alternative to doing a gaussian shape retrieval
  (e.g., see StateMappingGaussian). We have fewer retrieved variables
  for the number of levels we run the forward model on.

  For additional information see docs for StateMapping class.
*******************************************************************/
class StateMappingBasisMatrix : public StateMapping  {
public:
  StateMappingBasisMatrix(const blitz::Array<double, 2> Basis_matrix);
  virtual ~StateMappingBasisMatrix() {}
  
  virtual ArrayAd<double, 1> mapped_state(const ArrayAd<double, 1>& retrieval_values) const;
  virtual ArrayAd<double, 1> retrieval_state
  (const ArrayAd<double, 1>& initial_values) const;
  virtual std::string name() const { return "basis matrix"; }
  virtual boost::shared_ptr<StateMapping> clone() const
  { return boost::make_shared<StateMappingBasisMatrix>(basis_matrix()); }

  //-----------------------------------------------------------------------
  /// Matrix that takes us from the retrieval grid to the forward
  /// model grid. This will be m x n, with m the size of the forward
  /// grid and n the size of the retrieval grid (so m > n).
  ///
  /// When comparing with py_retrieve, note that this is the transpose
  /// of the basis matrix it uses.
  //-----------------------------------------------------------------------
    
  const blitz::Array<double, 2>& basis_matrix() const
  { return basis_matrix_; }
  
  //-----------------------------------------------------------------------
  /// Inverse of basis matrix, which goes from the forward model grid
  /// to the model grid.
  //-----------------------------------------------------------------------
  const blitz::Array<double, 2>& inverse_basis_matrix() const
  { return inverse_basis_matrix_;}
private:
  blitz::Array<double, 2> basis_matrix_, inverse_basis_matrix_;
  friend class boost::serialization::access;
  StateMappingBasisMatrix() {}
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateMappingBasisMatrix);

#endif
