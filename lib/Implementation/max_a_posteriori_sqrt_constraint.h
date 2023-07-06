#ifndef MAX_A_POSTERIORI_SQRT_CONSTRAINT_H
#define MAX_A_POSTERIORI_SQRT_CONSTRAINT_H
#include <max_a_posteriori.h>
#include <model_measure_standard.h>

namespace FullPhysics {
/******************************************************************
  This is a variation of MaxAPosterioriStandard, where we supply
  a "sqrt_constraint" to use instead of the a_priori matrix.

  Our original MaxAPosterioriStandard uses a Cholesky decomposition,
  while py-retrieve uses a different decomposition based on SVD
  (there are numerous decomposition of B = A A^T).
  
  The two are in some sense equivalent - the norm of the residual
  for the apriori augmentation is the same. However the specific
  values of the residual are different - basically these are rotations
  of each other.

  I'm not sure that it matters, but there may be calculation in
  py-retrieve that depends on the specific order of the residual
  relative to its sqrt_constraint matrix. So we have this variation
  that uses a supplied sqrt_constraint matrix rather than taking an 
  a priori matrix. Note that the sqrt_constraint we take is the
  transpose of what py-retrieve has - for some reason it transposes
  this after calculating it.

  Note by the way that despite the name, sqrt_constraint *isn't*
  actually the sqrt of the constraint matrix. The

  constraint_matrix = a_priori_cov^-1

  And the matrix is decomposed so

  sqrt_constraint * sqrt_constraint^T = constraint_matrix

  The decomposition by something like scipy.linalg.sqrtm is also
  a common way to handle the a priori augmentation, but this has

  r * r = constraint_matrix

  if r = sqrtm(constraint_matrix).

  sqrt_constraint is *not* equal to r, py-retrieve doesn't actually
  do a sqrtm.

  Instead the sqrt_constraint is like a cholesky decomposition of the
  constraint_matrix, it just isn't confined to being a lower
  triangular matrix. I don't believe the decomposition is unique
  either (although I'm not sure about that) - but since we take the
  decomposition as an input there isn't any ambiguity in this class.

  Note if you use this class then the calls to a_priori_cov_chol_inv()
  and a_priori_cov_chol() actually returns sqrt_constraint and
  sqrt_constraint^-1. I briefly considered overriding these functions
  to raise an Exception, since these aren't actually lower triangular
  matrices you get from a Cholesky decomposition. But I decide against
  that, for most uses I could think of you probably don't actually
  care that this is a different decomposition. But we can reconsider
  that if needed in the future - we really do have function names here
  that lie.
*******************************************************************/
class MaxAPosterioriSqrtConstraint : 
  virtual public MaxAPosteriori, virtual public ModelMeasureStandard {
public:
  MaxAPosterioriSqrtConstraint(const boost::shared_ptr<ForwardModel>& fm,
          const boost::shared_ptr<Observation>& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> sqrt_constraint);

  MaxAPosterioriSqrtConstraint(const std::vector<boost::shared_ptr<ForwardModel> >& fm,
          const std::vector<boost::shared_ptr<Observation> >& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> sqrt_constraint);
  
  virtual ~MaxAPosterioriSqrtConstraint() {}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "MaxAPosterioriSqrtConstraint"; }

private:
  MaxAPosterioriSqrtConstraint() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(MaxAPosterioriSqrtConstraint);
#endif
