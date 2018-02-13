#include <nlls_max_a_posteriori.h>


using namespace FullPhysics;
using namespace blitz;



#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NLLSMaxAPosteriori, CostFunc)
.def(luabind::constructor< const boost::shared_ptr<MaxAPosteriori>& >())
.def("parameters", (void( NLLSMaxAPosteriori::*)(const blitz::Array<double, 1>&))&NLLSMaxAPosteriori::parameters)
.def("max_a_posteriori", &NLLSMaxAPosteriori::max_a_posteriori)
REGISTER_LUA_END()
#endif



Array<double, 1> NLLSMaxAPosteriori::residual()
{
  if(R.size() <= 0) {

    assert_parameter_set_correctly();
//    MAP->parameters(X);  // If everything implemented correctly,
                                // there should be no need for this.

    //  For some mathematical models it is just practical to compute
    //  the model and its Jacobian at the same time.  Here, we want to
    //  check whether or not the model assigned to this NLLS adaptor
    //  (MAPE-problem to NLLS-problem) is implemented such that the
    //  model and its Jacobian are computed simultaneously.
    //
    bool j_computed_before = MAP->jacobean_computed();
    MAP->model_eval();
    bool j_computed_after = MAP->jacobean_computed();

    //  If the model Jacobian is computed as well when the model
    //  is computed, here we may as well compute the NLLS problem 
    //  residual Jacobian, which is not the same as the model Jacobian.  
    //  Usually the computation of the model Jacobian is much more
    //  expensive than the computation of the residual Jacobian based
    //  on the model Jacobian.
    //
    if( (!j_computed_before) && j_computed_after ) {
      increment_num_der1_evaluations();
      J.reference(MAP->weighted_jacobian_aug());
    }

    //  Compute the NLLS problem residual.
    //
    increment_num_cost_evaluations();  
    R.reference(MAP->weighted_model_measure_diff_aug());

  }
  return R.copy();
}


Array<double, 2> NLLSMaxAPosteriori::jacobian()
{
  if(J.size() <= 0) {

    assert_parameter_set_correctly();
//    MAP->parameters(X);

    //  For some mathematical models it is just practical to compute
    //  the model and its Jacobian at the same time.  Here, we want to
    //  check whether or not the model assigned to this NLLS adaptor
    //  (MAPE-problem to NLLS-problem) is implemented such that the
    //  model and it Jacobian are computed simultaneously.
    //
    bool m_computed_before = MAP->model_computed();
    MAP->jacobian_eval();
    bool m_computed_after = MAP->model_computed();

    //  If the model is computed as well when the model Jacobian
    //  is computed, here we may as well compute the NLLS problem 
    //  residual, which is not the same as the model.  Usually the
    //  computation of the model is much more expensive than the
    //  computation of the residual based on the model.
    //
    if( (!m_computed_before) && m_computed_after ) {
      increment_num_cost_evaluations();
      R.reference(MAP->weighted_model_measure_diff_aug());
    }

    //  Compute the NLLS problem residual Jacobian.
    //
    increment_num_der1_evaluations();
    J.reference(MAP->weighted_jacobian_aug());

  }
  return J.copy();
}


void NLLSMaxAPosteriori::parameters(const blitz::Array<double, 1>& x)
{
  MAP->parameters(x);
  NLLSProblemState::parameters(x);
}
