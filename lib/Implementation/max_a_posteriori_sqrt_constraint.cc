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
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MaxAPosteriori)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableMaxAPosterioriSqrtConstraint)
    & FP_NVP(K_x) & FP_NVP(msrmnt_jacobian_x) & FP_NVP_(mapping);
}

FP_IMPLEMENT(MaxAPosterioriSqrtConstraint);
FP_OBSERVER_SERIALIZE(MaxAPosterioriSqrtConstraint);
#endif

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

MaxAPosterioriSqrtConstraint::MaxAPosterioriSqrtConstraint
(const boost::shared_ptr<ForwardModel>& forward_model,
 const boost::shared_ptr<Observation>& observation, 
 const boost::shared_ptr<StateVector>& state_vector,
 const Array<double, 1> a_priori_params,
 const blitz::Array<double, 2> sqrt_constraint,
 const boost::shared_ptr<StateMapping>& in_map)
: ModelMeasureStandard(forward_model, observation, state_vector),
  mapping_(in_map)
{
  blitz::firstIndex i1; blitz::secondIndex i2; blitz::thirdIndex i3;
  
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
  blitz::Array<double, 1> x(sv->observer_claimed_size());
  x(blitz::Range::all()) = 0;
  blitz::Array<double, 1> z = mapping_->retrieval_state(x).value();
  if(Xa.rows() != z.rows())
    throw Exception("A priori state vector size and retrieval vector size expected by the model are not equal. :( ");
}

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

MaxAPosterioriSqrtConstraint::MaxAPosterioriSqrtConstraint
(const std::vector<boost::shared_ptr<ForwardModel> >& forward_model,
 const std::vector<boost::shared_ptr<Observation> >& observation, 
 const boost::shared_ptr<StateVector>& state_vector,
 const Array<double, 1> a_priori_params,
 const blitz::Array<double, 2> sqrt_constraint,
 const boost::shared_ptr<StateMapping>& in_map)
: ModelMeasureStandard(forward_model, observation, state_vector),
  mapping_(in_map)
{
  blitz::firstIndex i1; blitz::secondIndex i2; blitz::thirdIndex i3;
  
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
  blitz::Array<double, 1> x(sv->observer_claimed_size());
  x(blitz::Range::all()) = 0;
  blitz::Array<double, 1> z = mapping_->retrieval_state(x).value();
  if(Xa.rows() != z.rows())
    throw Exception("A priori state vector size and retrieval vector size expected by the model are not equal. :( ");
}

void MaxAPosterioriSqrtConstraint::parameters(const blitz::Array<double, 1>& z)
{
  if(!parameters_different(z)) return;
  ModelMeasure::parameters(z);
  blitz::Array<double, 1> x = mapping_->mapped_state(z).value();
  sv->update_state(x);
  notify_update_do(*this);
}

blitz::Array<double, 2> MaxAPosterioriSqrtConstraint::model_measure_diff_jacobian_fm()
{
  Array<double, 2> mjac(measurement_jacobian_fm());
  Array<double, 2> jac(jacobian_fm());
  Array<double, 2> result(measurement_size(), jac.cols());
  if(!msrmnt_is_const) {
    if(jac.rows() != mjac.rows() ||
       jac.cols() != mjac.cols())
      throw Exception("Model jacobian_x and measurement jacobian_x need to be the same size");
    result = jac - mjac;
  } else {
    result = jac;
  }
  return result;
}

void MaxAPosterioriSqrtConstraint::measurement_eval()
{
  if(msrmnt.size() > 0)
    return;
  blitz::Array<double, 1> z = parameters();
  ModelMeasureStandard::measurement_eval();
  if(!msrmnt_is_const) {
    msrmnt_jacobian_x.reference(msrmnt_jacobian.copy());
    msrmnt_jacobian.reference(mapping_->jacobian_retrieval(z, msrmnt_jacobian_x));
  }
}

void MaxAPosterioriSqrtConstraint::radiance_from_fm(bool skip_check)
{
  // ModelMeasureStandard gets the model and jacobian value for the
  // full state grid. Go ahead and call that, skipping the check on M
  // and K since they aren't the right size here.
  blitz::Array<double, 1> z = parameters();
  ModelMeasureStandard::radiance_from_fm(true);
  K_x.reference(K.copy());
  K.reference(mapping_->jacobian_retrieval(z, K_x));
  if(!skip_check)
    assert_model_correct(M);
  if(!skip_check)
    assert_jacobian_correct(K);
}

