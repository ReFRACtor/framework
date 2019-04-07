#include <meyer_nlls_problem.h>
#include <unit_test_support.h>
#include <fp_exception.h>
#include <nlls_solver_gsl_sm.h>


using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(nlls_solver_gsl_sm_meyer, GlobalFixture)

/* convergence check thresholds */
double x_tol=1.0e-6, g_tol=6.0555e-06, f_tol=1.0e-6;
bool verbose=false;


BOOST_AUTO_TEST_CASE(meyer_gsl_sm__default)
{
  //  Create the problem
  //
  Array<double, 1> x0(3); x0 = 0.02, 4000.0, 250.0;
  boost::shared_ptr<MeyerNLLSProblem> nlls(new MeyerNLLSProblem);

  //  Set the problem at the initial guess.
  //
  nlls->parameters(x0);

  //  Create the solver with almost all default setting,
  //  and solve the NLLS problem.
  //
  NLLSSolverGSLSM solver(nlls, 200);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)NLLSSolverGSLSM::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

  // std::cout 
  //    << "Testing NLLSSolverGSLSM with Meyer function:" << std::endl
  //    << "   Number of residual function evaluations = " << n_f_calls << std::endl
  //    << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
  //    << "   Final solver status = " << solver.status_str() << std::endl
  //    << "   Final problem status (point) = " << nlls->parameters() << std::endl
  //    << "   Final problem status (cost value) = " << cst << std::endl
  //    << "   Final problem status (gradient) = " << nlls->gradient() << std::endl;
  // for( int i=0; i<=solver.num_accepted_steps(); i++ )
  //    std::cout 
  //       << "   ========================================" << std::endl
  //       << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
  //       << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl
  //       << "   Grad["<<i<<"] = " << solver.gradient_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.gradient_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);
  BOOST_CHECK(fabs(sum(abs(nlls->gradient()-solver.gradient_at_accepted_points()[iLast]))) < 1e-8);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)NLLSSolverGSLSM::SUCCESS);
  BOOST_CHECK(n_f_calls < 155);
  BOOST_CHECK(n_j_calls <= n_f_calls);
  BOOST_CHECK_CLOSE(cst, 43.9729, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 0.0056096, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(1), 6181.35, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(2), 345.224, 0.01);
}


BOOST_AUTO_TEST_CASE(meyer_gsl_sm__lm_more_svd)
{
  //  Create the problem
  //
  Array<double, 1> x0(3); x0 = 0.02, 4000.0, 250.0;
  boost::shared_ptr<MeyerNLLSProblem> nlls(new MeyerNLLSProblem);

  //  Set the problem at the initial guess.
  //
  nlls->parameters(x0);

  //  Create the solver and solve the NLLS problem.
  //  Solver is set to the Lev/Mar trust region method
  //  with More's scaling, and it will use SVD for solving 
  //  linear systems.
  //
  gsl_multifit_nlinear_parameters fdf_params=gsl_multifit_nlinear_default_parameters();
  fdf_params.trs = gsl_multifit_nlinear_trs_lm;
  fdf_params.scale = gsl_multifit_nlinear_scale_more;
  fdf_params.solver = gsl_multifit_nlinear_solver_svd;
  //
  NLLSSolverGSLSM solver(nlls, 200, fdf_params);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)NLLSSolverGSLSM::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

  // std::cout 
  //    << "Testing NLLSSolverGSLSM with Meyer function:" << std::endl
  //    << "   Number of residual function evaluations = " << n_f_calls << std::endl
  //    << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
  //    << "   Final solver status = " << solver.status_str() << std::endl
  //    << "   Final problem status (point) = " << nlls->parameters() << std::endl
  //    << "   Final problem status (cost value) = " << cst << std::endl
  //    << "   Final problem status (gradient) = " << nlls->gradient() << std::endl;
  // for( int i=0; i<=solver.num_accepted_steps(); i++ )
  //    std::cout 
  //       << "   ========================================" << std::endl
  //       << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
  //       << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl
  //       << "   Grad["<<i<<"] = " << solver.gradient_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.gradient_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);
  BOOST_CHECK(fabs(sum(abs(nlls->gradient()-solver.gradient_at_accepted_points()[iLast]))) < 1e-8);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)NLLSSolverGSLSM::SUCCESS);
  BOOST_CHECK(n_f_calls < 160);
  BOOST_CHECK(n_j_calls <= n_f_calls);
  BOOST_CHECK_CLOSE(cst, 43.9729, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 0.0056096, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(1), 6181.35, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(2), 345.224, 0.01);
}


BOOST_AUTO_TEST_CASE(meyer_gsl_sm__lmaccel_more_svd)
{
  //  Create the problem
  //
  Array<double, 1> x0(3); x0 = 0.02, 4000.0, 250.0;
  boost::shared_ptr<MeyerNLLSProblem> nlls(new MeyerNLLSProblem);

  //  Set the problem at the initial guess.
  //
  nlls->parameters(x0);

  //  Create the solver and solve the NLLS problem.  Solver
  //  is set to the Lev/Mar accelerated trust region method
  //  with More's scaling, and it will use SVD for solving 
  //  linear systems.
  //
  gsl_multifit_nlinear_parameters fdf_params=gsl_multifit_nlinear_default_parameters();
  fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
  fdf_params.scale = gsl_multifit_nlinear_scale_more;
  fdf_params.solver = gsl_multifit_nlinear_solver_svd;
  fdf_params.fdtype =  GSL_MULTIFIT_NLINEAR_CTRDIFF;
  //
  NLLSSolverGSLSM solver(nlls, 200, fdf_params);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)NLLSSolverGSLSM::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

  // std::cout 
  //    << "Testing NLLSSolverGSLSM with Meyer function:" << std::endl
  //    << "   Number of residual function evaluations = " << n_f_calls << std::endl
  //    << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
  //    << "   Final solver status = " << solver.status_str() << std::endl
  //    << "   Final problem status (point) = " << nlls->parameters() << std::endl
  //    << "   Final problem status (cost value) = " << cst << std::endl
  //    << "   Final problem status (gradient) = " << nlls->gradient() << std::endl;
  // for( int i=0; i<=solver.num_accepted_steps(); i++ )
  //    std::cout 
  //       << "   ========================================" << std::endl
  //       << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
  //       << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl
  //       << "   Grad["<<i<<"] = " << solver.gradient_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.gradient_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);
  BOOST_CHECK(fabs(sum(abs(nlls->gradient()-solver.gradient_at_accepted_points()[iLast]))) < 1e-8);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)NLLSSolverGSLSM::SUCCESS);
  BOOST_CHECK(n_f_calls < 100);
  BOOST_CHECK(n_j_calls <= n_f_calls);
  BOOST_CHECK_CLOSE(cst, 43.9729, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 0.0056096, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(1), 6181.35, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(2), 345.224, 0.01);
}


BOOST_AUTO_TEST_CASE(meyer_gsl_sm__subspace2D_more_svd)
{
  //  Create the problem
  //
  Array<double, 1> x0(3); x0 = 0.02, 4000.0, 250.0;
  boost::shared_ptr<MeyerNLLSProblem> nlls(new MeyerNLLSProblem);

  //  Set the problem at the initial guess.
  //
  nlls->parameters(x0);

  //  Create the solver and solve the NLLS problem.
  //  Solver is set to the subspace2D trust region method
  //  with More's scaling, and it will use SVD for solving 
  //  linear systems.
  //
  gsl_multifit_nlinear_parameters fdf_params=gsl_multifit_nlinear_default_parameters();
  fdf_params.trs = gsl_multifit_nlinear_trs_subspace2D;
  fdf_params.scale = gsl_multifit_nlinear_scale_more;
  fdf_params.solver = gsl_multifit_nlinear_solver_svd;
  //
  NLLSSolverGSLSM solver(nlls, 200, fdf_params);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)NLLSSolverGSLSM::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

  // std::cout 
  //    << "Testing NLLSSolverGSLSM with Meyer function:" << std::endl
  //    << "   Number of residual function evaluations = " << n_f_calls << std::endl
  //    << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
  //    << "   Final solver status = " << solver.status_str() << std::endl
  //    << "   Final problem status (point) = " << nlls->parameters() << std::endl
  //    << "   Final problem status (cost value) = " << cst << std::endl
  //    << "   Final problem status (gradient) = " << nlls->gradient() << std::endl;
  // for( int i=0; i<=solver.num_accepted_steps(); i++ )
  //    std::cout 
  //       << "   ========================================" << std::endl
  //       << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
  //       << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl
  //       << "   Grad["<<i<<"] = " << solver.gradient_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.gradient_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);
  BOOST_CHECK(fabs(sum(abs(nlls->gradient()-solver.gradient_at_accepted_points()[iLast]))) < 1e-8);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)NLLSSolverGSLSM::SUCCESS);
  BOOST_CHECK(n_f_calls < 15);
  BOOST_CHECK(n_j_calls <= n_f_calls);
  BOOST_CHECK_CLOSE(cst, 43.9729, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 0.0056096, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(1), 6181.35, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(2), 345.224, 0.01);
}


BOOST_AUTO_TEST_CASE(meyer_gsl_sm__ddogleg_more_svd)
{
  //  Create the problem
  //
  Array<double, 1> x0(3); x0 = 0.02, 4000.0, 250.0;
  boost::shared_ptr<MeyerNLLSProblem> nlls(new MeyerNLLSProblem);

  //  Set the problem at the initial guess.
  //
  nlls->parameters(x0);

  //  Create the solver and solve the NLLS problem.
  //  Solver is set to the ddogleg trust region method
  //  with More's scaling, and it will use SVD for solving 
  //  linear systems.
  //
  gsl_multifit_nlinear_parameters fdf_params=gsl_multifit_nlinear_default_parameters();
  fdf_params.trs = gsl_multifit_nlinear_trs_ddogleg;
  fdf_params.scale = gsl_multifit_nlinear_scale_more;
  fdf_params.solver = gsl_multifit_nlinear_solver_svd;
  //
  NLLSSolverGSLSM solver(nlls, 200, fdf_params);
  BOOST_CHECK_EQUAL((int)solver.status(), (int)NLLSSolverGSLSM::UNTRIED);
  solver.solve();

  int n_f_calls = nlls->num_residual_evaluations();
  int n_j_calls = nlls->num_jacobian_evaluations();
  double cst = nlls->cost();

  // std::cout 
  //    << "Testing NLLSSolverGSLSM with Meyer function:" << std::endl
  //    << "   Number of residual function evaluations = " << n_f_calls << std::endl
  //    << "   Number of jacobian function evaluations = " << n_j_calls << std::endl
  //    << "   Final solver status = " << solver.status_str() << std::endl
  //    << "   Final problem status (point) = " << nlls->parameters() << std::endl
  //    << "   Final problem status (cost value) = " << cst << std::endl
  //    << "   Final problem status (gradient) = " << nlls->gradient() << std::endl;
  // for( int i=0; i<=solver.num_accepted_steps(); i++ )
  //    std::cout 
  //       << "   ========================================" << std::endl
  //       << "   At point["<<i<<"] " << solver.accepted_points()[i] << std::endl
  //       << "   Cost["<<i<<"] = " << solver.cost_at_accepted_points()[i] << std::endl
  //       << "   Grad["<<i<<"] = " << solver.gradient_at_accepted_points()[i] << std::endl;

  BOOST_CHECK_EQUAL((int) solver.accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.cost_at_accepted_points().size(), solver.num_accepted_steps()+1);
  BOOST_CHECK_EQUAL((int) solver.gradient_at_accepted_points().size(), solver.num_accepted_steps()+1);

  int iLast = solver.num_accepted_steps();
  BOOST_CHECK_EQUAL(cst, solver.cost_at_accepted_points()[iLast]);
  BOOST_CHECK_CLOSE(sum(abs(nlls->parameters()-solver.accepted_points()[iLast])), 0.0, 1e-12);
  BOOST_CHECK(fabs(sum(abs(nlls->gradient()-solver.gradient_at_accepted_points()[iLast]))) < 1e-8);

  BOOST_CHECK_EQUAL((int)solver.status(), (int)NLLSSolverGSLSM::SUCCESS);
  BOOST_CHECK(n_f_calls < 15);
  BOOST_CHECK(n_j_calls <= n_f_calls);
  BOOST_CHECK_CLOSE(cst, 43.9729, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(0), 0.0056096, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(1), 6181.35, 0.01);
  BOOST_CHECK_CLOSE(nlls->parameters()(2), 345.224, 0.01);
}



BOOST_AUTO_TEST_SUITE_END()
