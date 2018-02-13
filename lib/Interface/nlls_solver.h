#ifndef NLLS_SOLVER_H
#define NLLS_SOLVER_H
#include <iterative_solver_der.h>
#include <nlls_problem.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for the solvers of the 
///        Nonlinear-Least-Squares Problem.
///
/// This is the base class for all Nonlinear-Least-Squares
/// solvers.
///
/// This class is associated with a problem (NLLSProblem)
/// because the problem interface is determined:
///   - provide a point in the parameter space
///   - evaluate the residual function (a vector function)
///     at the point
///   - evaluate the Jacobian (a matrix function) of the
///     residual at the point
//-----------------------------------------------------------------------

class NLLSSolver : 
    public IterativeSolverDer {

public:


//-----------------------------------------------------------------------
/// \brief Constructor
/// 
/// \param[in] max_cost_function_calls 
///            read related base class comments
///
/// \param[in] p
///            The Nonlinear Least Squares problem
///
/// \param[in] vrbs
///            read related base class comments
//-----------------------------------------------------------------------

  NLLSSolver(const boost::shared_ptr<NLLSProblem>& p,
             int max_cost_function_calls, bool vrbs)
    : IterativeSolverDer(max_cost_function_calls, vrbs),
      P(p)
  {}


  virtual ~NLLSSolver() {}


//-----------------------------------------------------------------------
/// Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "NLLSSolver"; }


protected:

  boost::shared_ptr<NLLSProblem> P;

};
}
#endif
