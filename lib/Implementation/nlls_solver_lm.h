#ifndef NLLS_SOLVER_LM_H
#define NLLS_SOLVER_LM_H

#include <Eigen/Dense>
#include <nlls_solver.h>


namespace FullPhysics {

/******************************************************************
  \brief Homemade NLLS solver based on Levenberg-Marquardt
         algorithm.
 
  This class is the implementation of J. J. More's
  version of the Levenberg-Marquardt NLLS solver.
*******************************************************************/
class NLLSSolverLM :
    public NLLSSolver {

public:

  /******************************************************************
   * \brief Provides NLLSSolverLM with very algorithm specific
   *        parameters.
   *
   * This class is used to collect and provide a solver of NLLSSolverLM
   * type with some of the parameters that is not commonly used by
   * other types of solvers.  The constructor sets each parameter to
   * a reasonable value; however, the trust region radius "tr_rad" is
   * problem dependent and could optionally be reset to a value more
   * appropriate for the problem to be solved.
  *******************************************************************/
  class Options
  {
  public:
    Options()
      : min_W(1.0e-9), tr_rad_tol(0.1), tr_rad(1000.0), cr_ratio_tol(0.0001)
    {}

    double min_W;  ///< A positive number and the smallest value for the
                   ///< diagonal entries of the diagonal weight matrix.

    double tr_rad_tol;  ///< A positive number and a tolerance for comparing to the
                        ///< trust region radius. A step within this tolerance
                        ///< distance from the radius is considered to be on the radius.

    double tr_rad;  ///< A positive number and the initial trust region radius.

    double cr_ratio_tol;  ///< A positive number and
                          ///< tolerance on the ratio of the actual to the predicted
                          ///< reduction of the value of the cost function after taking
                          ///< a step.  A ratio greater than this tolerance results in
                          ///< accepting the step.  Otherwise, the step is rejected.
  };


//-----------------------------------------------------------------------
///  \brief Constructor
/// 
///  \param[in] p
///             The Nonlinear Least Squares problem
///  
///  \param[in] max_cost_function_calls 
///             read comments on NLLSSolver::NLLSSolver(int32_t,bool)
/// 
///  \param[in] opt
///             A collection of some of the solver controlling
///             parameters. See the comments on NLLSSolverLM::Options.
/// 
///  \param[in] dx_tol_abs
///             An absolute tolerance on step size.
/// 
///  \param[in] dx_tol_rel
///             A relative tolerance on step size.
/// 
///  \param[in] g_tol_abs
///             An absolute tolerance on gradient size.
/// 
///  \param[in] g_tol_rel
///             A relative tolerance on gradient size.
/// 
///  \param[in] vrbs
///             read comments on NLLSSolver::NLLSSolver(int32_t,bool)
//-----------------------------------------------------------------------
  NLLSSolverLM( const boost::shared_ptr<NLLSProblem>& p, int max_cost_function_calls, 
                const NLLSSolverLM::Options& opt=NLLSSolverLM::Options(),
                double dx_tol_abs=0.000001, double dx_tol_rel=0.000001,
                double g_tol_abs=0.000001, double g_tol_rel=6.0555e-06,
                bool vrbs=false );


//-----------------------------------------------------------------------
///  \brief Destructor
//-----------------------------------------------------------------------

  virtual ~NLLSSolverLM()
  {}


//-----------------------------------------------------------------------
///  \brief The method that solves the optimization problem.
/// 
///  Read the comments on IterativeSolver::solve().
//-----------------------------------------------------------------------

  virtual void solve();


//-----------------------------------------------------------------------
/// 
///  \brief Prints solver state.
/// 
///  The method prints some information about the current
///  state of the solver and the problem, for example,
///  the current point and the cost function value and its
///  gradient at the current point.
//-----------------------------------------------------------------------

  virtual void print_state(std::ostream &ostr=std::cout);


//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "NLLSSolverLM"; }


protected:

  virtual status_t iterate();

  static status_t test_dx_rel( const Eigen::VectorXd& dx, const Eigen::VectorXd& x, double dx_tol_rel );

  static status_t test_dx_abs( const Eigen::VectorXd& dx, double dx_tol_abs );

  static status_t test_dx( const Eigen::VectorXd& dx, const Eigen::VectorXd& x, double dx_tol_rel, double dx_tol_abs );

  static status_t test_grad_rel( const Eigen::VectorXd& g, const Eigen::VectorXd& x, double cost, double g_tol_rel );

  static status_t test_grad_abs( const Eigen::VectorXd& g, double g_tol_abs );

  double Dx_tol_abs;
  double Dx_tol_rel;
  double G_tol_abs;
  double G_tol_rel;

  Options Opt;

  double CR_ratio;
  double Lambda;

  Eigen::VectorXd W;

  Eigen::VectorXd Dx;

};  // class NLLSProblemLM

}  // namespace FullPhysics

#endif
