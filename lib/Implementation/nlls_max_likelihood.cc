#include <nlls_max_likelihood.h>
#include "fp_serialize_support.h"


using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void NLLSMaxLikelihood::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NLLSProblem)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NLLSProblemState)
    & FP_NVP(ML);
}

FP_IMPLEMENT(NLLSMaxLikelihood);
#endif


#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NLLSMaxLikelihood, CostFunc)
.def(luabind::constructor< const boost::shared_ptr<MaxLikelihood>& >())
.def("parameters", (void( NLLSMaxLikelihood::*)(const blitz::Array<double, 1>&))&NLLSMaxLikelihood::parameters)
.def("max_likelihood", &NLLSMaxLikelihood::max_likelihood)
REGISTER_LUA_END()
#endif



Array<double, 1> NLLSMaxLikelihood::residual()
{
  if(R.size() <= 0) { 

    assert_parameter_set_correctly();
//    ML.set_parameters(X);  // If everything implemented correctly, there should be no need for this.

    //  For some mathematical models it is just practical to compute
    //  the model and its Jacobian at the same time.  Here, we want to
    //  check whether or not the model assigned to this NLLS adaptor
    //  (MLE-problem to NLLS-problem) is implemented such that the
    //  model and its Jacobian are computed simultaneously.
    //
    bool j_computed_before = ML->jacobean_computed();
    ML->model_eval();
    bool j_computed_after = ML->jacobean_computed();

    //  If the model Jacobian is computed as well when the model
    //  is computed, here we may as well compute the NLLS problem 
    //  residual Jacobian, which is not the same as the model Jacobian.  
    //  Usually the computation of the model Jacobian is much more
    //  expensive than the computation of the residual Jacobian based
    //  on the model Jacobian.
    //
    if( (!j_computed_before) && j_computed_after ) {
      increment_num_der1_evaluations();
      J.reference(ML->uncert_weighted_jacobian());
    }

    //  Compute the NLLS problem residual.
    //
    increment_num_cost_evaluations();  
    R.reference(ML->uncert_weighted_model_measure_diff());

  }
  return R.copy();
}


Array<double, 2> NLLSMaxLikelihood::jacobian()
{
  if(J.size() <= 0) {

    assert_parameter_set_correctly();
//    ML.parameters(X);

    //  For some mathematical models it is just practical to compute
    //  the model and its Jacobian at the same time.  Here, we want to
    //  check whether or not the model assigned to this NLLS adaptor
    //  (MLE-problem to NLLS-problem) is implemented such that the
    //  model and it Jacobian are computed simultaneously.
    //
    bool m_computed_before = ML->model_computed();
    ML->jacobian_eval();
    bool m_computed_after = ML->model_computed();

    //  If the model is computed as well when the model Jacobian
    //  is computed, here we may as well compute the NLLS problem 
    //  residual, which is not the same as the model.  Usually the
    //  computation of the model is much more expensive than the
    //  computation of the residual based on the model.
    //
    if( (!m_computed_before) && m_computed_after ) {
      increment_num_cost_evaluations();
      R.reference(ML->uncert_weighted_model_measure_diff());
    }

    //  Compute the NLLS problem residual Jacobian.
    //
    increment_num_der1_evaluations();
    J.reference(ML->uncert_weighted_jacobian());

  }
  return J.copy();
}


void NLLSMaxLikelihood::parameters(const blitz::Array<double, 1>& x)
{
  ML->parameters(x);
  NLLSProblemState::parameters(x);
}
