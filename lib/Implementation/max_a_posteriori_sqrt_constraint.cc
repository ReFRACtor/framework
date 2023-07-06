#include "max_a_posteriori_sqrt_constraint.h"
#include "linear_algebra.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void MaxAPosterioriSqrtConstraint::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ModelMeasureStandard)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MaxAPosteriori);
}

FP_IMPLEMENT(MaxAPosterioriSqrtConstraint);
#endif

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

MaxAPosterioriSqrtConstraint::MaxAPosterioriSqrtConstraint
(const boost::shared_ptr<ForwardModel>& forward_model,
 const boost::shared_ptr<Observation>& observation, 
 const boost::shared_ptr<StateVector>& state_vector,
 const Array<double, 1> a_priori_params,
 const blitz::Array<double, 2> sqrt_constraint)
: ModelMeasureStandard(forward_model, observation, state_vector) 
{
  blitz::firstIndex i1; blitz::secondIndex i2; blitz:thirdIndex i3;
  
  if(a_priori_params.rows() <= 0)
    throw Exception("The size of the a priori parameter values array is zero.");
  if(sqrt_constraint.rows() != sqrt_constraint.cols() )
    throw Exception("The a sqrt_constraint matrix must be a square matrix.");
  if(sqrt_constraint.rows() != a_priori_params.rows() )
    throw Exception("The number of rows and columns of sqrt_constraint matrix must equal the size of the a priori parameter values array.");
  
  Xa.reference(a_priori_params.copy());
  // Not really the cholesky inverse, but other than the name acts
  // the same in MaxAPosteriori.
  Sa_chol_inv.reference(sqrt_constraint.copy());
  Sa_chol.reference(generalized_inverse(Sa_chol_inv, 1e-20));
  Sa.resize(Xa.rows(), Xa.rows());
  Sa = blitz::sum(Sa_chol(i1, i3) * Sa_chol(i2,i3), i3);
  if(Xa.rows() != sv->observer_claimed_size())
    throw Exception("A priori state vector size and state vector size expected by the model are not equal. :( ");
}

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

MaxAPosterioriSqrtConstraint::MaxAPosterioriSqrtConstraint
(const std::vector<boost::shared_ptr<ForwardModel> >& forward_model,
 const std::vector<boost::shared_ptr<Observation> >& observation, 
 const boost::shared_ptr<StateVector>& state_vector,
 const Array<double, 1> a_priori_params,
 const blitz::Array<double, 2> sqrt_constraint)
: ModelMeasureStandard(forward_model, observation, state_vector) 
{
  blitz::firstIndex i1; blitz::secondIndex i2; blitz:thirdIndex i3;
  
  if(a_priori_params.rows() <= 0)
    throw Exception("The size of the a priori parameter values array is zero.");
  if(sqrt_constraint.rows() != sqrt_constraint.cols() )
    throw Exception("The a sqrt_constraint matrix must be a square matrix.");
  if(sqrt_constraint.rows() != a_priori_params.rows() )
    throw Exception("The number of rows and columns of sqrt_constraint matrix must equal the size of the a priori parameter values array.");
  
  Xa.reference(a_priori_params.copy());
  // Not really the cholesky inverse, but other than the name acts
  // the same in MaxAPosteriori.
  Sa_chol_inv.reference(sqrt_constraint.copy());
  Sa_chol.reference(generalized_inverse(Sa_chol_inv, 1e-20));
  Sa.resize(Xa.rows(), Xa.rows());
  Sa = blitz::sum(Sa_chol(i1, i3) * Sa_chol(i2,i3), i3);
  if(Xa.rows() != sv->observer_claimed_size())
    throw Exception("A priori state vector size and state vector size expected by the model are not equal. :( ");
}
