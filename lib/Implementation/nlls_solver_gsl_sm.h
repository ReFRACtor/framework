#ifndef NLLS_SOLVER_GSL_SM_H
#define NLLS_SOLVER_GSL_SM_H

#include <gsl/gsl_multifit_nlinear.h>
#include <nlls_solver.h>


namespace FullPhysics {


/******************************************************************
  This class is the base class for the solvers of the NLLS
  problem that are small-to-medium is size.  It requires the 
  newer GSL library that has solvers for small-to-medium and
  large problem.
*******************************************************************/
class NLLSSolverGSLSM :
    public NLLSSolver {

public:

//-----------------------------------------------------------------------
/// Initializes the solver.
///
/// \param[in] p
///            The Nonlinear Least Squares problem
/// 
/// \param[in] max_cost_function_calls 
///            read comments on NLLSSolver::NLLSSolver(int32_t,bool)
///
/// \param[in] fdf_params
///            .
///
/// \param[in] x_tol
///            Used in testing for a small step size relative to the 
///            current parameter vector.  The valuse should be 1.0^-d,
///            where d is the desired number of the accurate digits
///            in the minimizer.
///
/// \param[in] g_tol
///            Used in testing for a small gradient.
///
/// \param[in] f_tol
///            .
///
/// \param[in] vrbs
///            read comments on NLLSSolver::NLLSSolver(int32_t,bool)
//-----------------------------------------------------------------------
  NLLSSolverGSLSM( const boost::shared_ptr<NLLSProblem>& p, int32_t max_cost_function_calls, 
                   gsl_multifit_nlinear_parameters fdf_params=gsl_multifit_nlinear_default_parameters(),
                   double x_tol=1.0e-6, double g_tol=6.0555e-06, double f_tol=0.0, bool vrbs=false )
    : NLLSSolver(p,max_cost_function_calls,vrbs),
      FDF_params(fdf_params),
      X_tol(x_tol), G_tol(g_tol), F_tol(f_tol)
  {}


  virtual ~NLLSSolverGSLSM()
  {}


  virtual void solve();


//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "NLLSSolverGSLSM"; }


protected:


  double X_tol;
  double G_tol;
  double F_tol;

  gsl_multifit_nlinear_parameters FDF_params;


  virtual const gsl_multifit_nlinear_type* get_gsl_multifit_nlinear_solver()
  { return gsl_multifit_nlinear_trust; /*default*/ }


};  // class NLLSProblemGSLSM

}


#endif
