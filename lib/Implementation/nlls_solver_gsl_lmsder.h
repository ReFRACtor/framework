#ifndef NLLS_SOLVER_GSL_LMSDER_H
#define NLLS_SOLVER_GSL_LMSDER_H
#include <nlls_solver_gsl.h>

namespace FullPhysics {
/******************************************************************
  This class is the implementation of J. J. More's
  version of the Levenberg-Marquardt NLLS solver.
*******************************************************************/
class NLLSSolverGSLLMSDER : 
    virtual public NLLSSolverGSL {

public:

//-----------------------------------------------------------------------
/// Initializes the solver.
/// 
/// \param p The problem Input value
/// \param max_cost_function_calls Input value
/// \param dx_tol_abs Input value
/// \param dx_tol_rel Input value
/// \param g_tol Input value
/// \param vrbs Input value
//-----------------------------------------------------------------------

  NLLSSolverGSLLMSDER(const boost::shared_ptr<NLLSProblem>& p, int max_cost_function_calls, 
                      double dx_tol_abs=0.000001, double dx_tol_rel=0.000001, double g_tol=6.0555e-06, 
                      bool vrbs=false)
    : NLLSSolverGSL(p, max_cost_function_calls, dx_tol_abs, dx_tol_rel, g_tol, vrbs)
  {}

  virtual ~NLLSSolverGSLLMSDER() {}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "NLLSSolverGSLLMSDER"; }

protected:

  virtual const gsl_multifit_fdfsolver_type* get_gsl_multifit_fdfsolver()
  { return gsl_multifit_fdfsolver_lmsder; /*same as the base class*/ }
  NLLSSolverGSLLMSDER() {}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(NLLSSolverGSLLMSDER);
#endif
