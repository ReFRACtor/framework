#ifndef NLLS_SOLVER_LM_H
#define NLLS_SOLVER_LM_H

#include <Eigen/Dense>
#include <nlls_solver.h>


namespace FullPhysics {

/**************************************************************//**
  \brief Homemade NLLS solver based on Levenberg-Marquardt
         algorithm.
 
  This class is the implementation of J. J. More's
  version of the Levenberg-Marquardt NLLS solver.
*******************************************************************/

class NLLSSolverLM :
    public NLLSSolver {

public:


//-----------------------------------------------------------------------
///  \brief Provides NLLSSolverLM with algorithm specific input parameters.
/// 
///  This class is used to collect and provide a solver of NLLSSolverLM
///  type with some of the parameters that is not commonly used by
///  other types of solvers.  The constructor sets each parameter to
///  a reasonable value; however, the trust region radius "tr_rad" is
///  problem dependent and could optionally be reset to a value more
///  appropriate for the problem to be solved.
//-----------------------------------------------------------------------

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
///  \brief For accumulating some NLLSSolverLM algorithm intermediate 
///         parameters.
/// 
///  The purpose of this class is to accumulate some solver intermediate
///  parameter values that could change from one iteration to the next.  
///  The values are accumulated either for recording purpose or for use 
///  by the solver for the next iteration.
///  
///  This parameters also give a measure of the linearity of the cost
///  function in the direction of taken steps.
//-----------------------------------------------------------------------

  class AlgParams
  {
  public:
    AlgParams()
      : cr_ratio(0.0), lambda(0.0), tr_rad(0.0), scaled_step_norm(0.0)
    {}

    double cr_ratio;  ///< The ratio of the actual reduction in the value of the
                      ///< cost function to the predicted reduction in the value of
                      ///< the cost function due to a step.  A value close to one
                      ///< suggest that the cost function is close to be linear
                      ///< over the step.  (needed from one iteration to the
                      ///< next and for recording purpose)

    double lambda;  ///< The Levenberg/Marquardt parameter: Its value is modified to 
                    ///< control the size and the direction of the step.  If the value
                    ///< of lambda is increased, then the step becomes smaller, and it 
                    ///< gets turned away from the Gauss-Newton step direction and 
                    ///< turned towards steepest descent direction.  (needed from one
                    ///< iteration to the next and for recording purpose)

    double tr_rad;  ///< This is the radius of a trust region centered at
                    ///< the current point where the algorithm determines
                    ///< that a linear approximation to the cost function
                    ///< is a reasonable approximation.  (needed from one 
                    ///< iteration to the next and for recording purpose)

    double scaled_step_norm;  ///< The norm of the scaled step taken at the current point.
                              ///< (for recording purpose only)

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
///  \brief Returns a std::vector<AlgParams> of intermediate algorithm
///         parameters.
/// 
///  This method returns a std vector of intermediate algorithm
///  parameters (AlgParams) recorded after each accepted step.  The
///  recording can be done more frequently, for example, once per
///  any attempted step (accepted or rejected); however, with the 
///  current implementation, the recording is once for each accepted
///  step.  Therefore, 
///    - algorithm_parameters().size() == IterativeSolver::num_accepted_steps()
/// 
///  With the current implementation, algorithm_parameters()[i] refers
///  to the collection of the intermediate parameters generated when
///  computing the accepted step that takes the problem from
///  accepted_points()[i] to accepted_points()[i+1] for
///  0 <= i < num_accepted_steps().
/// 
///  \return A vector of AlgParams
//-----------------------------------------------------------------------

  virtual const std::vector<AlgParams>& algorithm_parameters() const
  { return Alg_params_for_each_step; }


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

  static status_t test_grad( const Eigen::VectorXd& g, const Eigen::VectorXd& x, double cost, double g_tol_rel, double g_tol_abs );


//-----------------------------------------------------------------------
///  \brief Prints some solver intermediate parameters.
/// 
///  The method prints some intermediate solver parameters that
///  help to understand how linear the cost function is locally.
///  This method can be public, but it is just not public to 
///  indicate the parameters are very specific to this algorithm.
//-----------------------------------------------------------------------

  void print_state_linearity(std::ostream &ostr=std::cout);


//-----------------------------------------------------------------------
///  \brief For recording intermediate algorithm parameters
/// 
///  This method is called to record some useful intermediate algorithm
///  parameters.  They can be studied to learn about the performance of
///  the Lev/Mar iterative solver.
/// 
///  This method can be called once per any attempted step (accepted or 
///  rejected).  However, in the current implementation of the code, this
///  method is called once after each accepted step to record the parametes
///  for accepted steps only.
/// 
///  \param[in] ap
///             NLLSProblemLM::AlgParams type object containing 
///             intermediate algorithm parameters
//-----------------------------------------------------------------------

  void record_alg_params_after_step(const AlgParams& ap)
  { Alg_params_for_each_step.push_back(ap); }


  double Dx_tol_abs;
  double Dx_tol_rel;
  double G_tol_abs;
  double G_tol_rel;

  const Options Opt;
  AlgParams Ap;
  Eigen::VectorXd W;
  Eigen::VectorXd Dx;


private:

  std::vector<AlgParams> Alg_params_for_each_step;  ///< For recording intermediate algorithm parameters
                                                    ///< throughout the iterative process.

};  // class NLLSProblemLM

}  // namespace FullPhysics

#endif
