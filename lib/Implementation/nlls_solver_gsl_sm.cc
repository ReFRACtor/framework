#include <gsl/gsl_blas.h>
#include <fp_gsl_matrix.h>
#include <gsl_sm_lsp.h>
#include <nlls_solver_gsl_sm.h>


using namespace FullPhysics;
using namespace blitz;


void print_state( uint32_t iter,  gsl_multifit_nlinear_workspace* s, int status )
{
  double c = gsl_blas_dnrm2( gsl_multifit_nlinear_residual(s) );
  printf( "Solver '%s/%s';  iter = %3u;  (|f(x)|^2)/2 = %g;  status = %s\n", 
          gsl_multifit_nlinear_name(s), gsl_multifit_nlinear_trs_name(s),
          iter, c*c/2.0, gsl_strerror(status) );

  printf( "Where x is\n" );
  (void) gsl_vector_fprintf(stdout, gsl_multifit_nlinear_position(s), "%25.14lf");

//  printf( "The gradient g(x) is\n");
//  (void) gsl_vector_fprintf(stdout, s->g, "%25.14lf");

  printf( "\n" );
}


void NLLSSolverGSLSM::solve()
{
  //  Selecting a problem
  gsl_multifit_nlinear_fdf f = gsl_sm_get_lsp_fdf(&(*P));

  //  select the solver
  const gsl_multifit_nlinear_type * T = get_gsl_multifit_nlinear_solver();

  //  for gsl status
  int gsl_status = GSL_FAILURE;
  int gsl_status_conv = GSL_CONTINUE;

  gsl_multifit_nlinear_workspace *s = gsl_multifit_nlinear_alloc(T, &FDF_params, f.n, f.p);

  blitz::Array<double, 1> X(P->parameters());

  uint32_t num_step = 0;
  stat = UNTRIED;
  if( s )
    if( !(gsl_status = gsl_multifit_nlinear_init(GslVector(X).gsl(), &f, s)) ) {
      stat = CONTINUE;
      do {
        num_step++;

        //  The following three lines are only for recording purpose.
        //  Info at the initial guess (the starting point) is also 
        //  recorded here.
        //
        record_cost_at_accepted_point(P->cost());
        record_accepted_point(P->parameters());
        record_gradient_at_accepted_point(P->gradient());

        //  A return status of GSL_ENOPROG by the following function
        //  does not suggest convergence, and it only means that the
        //  function call did not encounter an error.  However, a
        //  return status of GSL_ENOPROG means that the function call
        //  did not make any progress.  Stalling can happen at a true
        //  minimum or some other reason.
        //
        gsl_status = gsl_multifit_nlinear_iterate(s);
        if( (gsl_status != GSL_SUCCESS) && (gsl_status != GSL_ENOPROG) ) {
          stat = ERROR;
          break;
        }

        //  The gradient check performed by the following GSL provided
        //  function is a scaled gradient check.  First the gradient 
        //  is scaled such that it is unitless, and then the convergence
        //  check is performed.  Scaling the gradient for convergence
        //  check helps to avoid some problems, but it results in other
        //  problem.
        //  
        int info=0;
        gsl_status_conv = gsl_multifit_nlinear_test(X_tol, G_tol, F_tol, &info, s);

        //  This code segment is a temporary fix.  The above GSL
        //  provided convergence test routine does not work correctly
        //  when the value of the cost function is very large.
        //  This code segment (if statement) is written such that it
        //  can be deleted without causing error or requiring any
        //  change to the rest of the code.  The effect of the code
        //  segment is that both the scaled gradient test (above) and
        //  the unscaled gradient test (below) should satisfy the
        //  gradient-based convergence tests in order to assume
        //  convergence based on gradient test.
        //
        if(info == 2) {
          if(sum(abs(P->gradient())) > G_tol) {
            gsl_status_conv = GSL_CONTINUE;
            info = 0;
          }
        }

        if( (info == 1) || (info == 2) ) {
          stat = SUCCESS;
          break;
        }

        if(gsl_status == GSL_ENOPROG) {
          stat = STALLED;
          break;
        }

        if( verbose ) print_state(num_step, s, gsl_status_conv);
      } while ( (gsl_status_conv == GSL_CONTINUE)
                && (P->num_cost_evaluations() < max_cost_f_calls)
                && (P->message() == NLLSProblem::NONE) );
    }

  //  The following three lines are only for recording purpose.
  //
  record_cost_at_accepted_point(P->cost());
  record_accepted_point(P->parameters());
  record_gradient_at_accepted_point(P->gradient());

  if( verbose && (stat != CONTINUE) )
    print_state( num_step, s, ((stat == ERROR)?gsl_status:gsl_status_conv) );

  if( s ) gsl_multifit_nlinear_free(s);
}
