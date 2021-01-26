#include <fp_gsl_matrix.h>
#include <fp_exception.h>
#include "fp_serialize_support.h"
#include <gsl_mdm.h>
#include <cost_minimizer_gsl.h>
#include<cmath>


using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void CostMinimizerGSL::serialize(Archive & ar,
				 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CostMinimizer)
    & FP_NVP(Size_tol) & FP_NVP(Initial_step_size);
}

FP_IMPLEMENT(CostMinimizerGSL);
#endif



boost::shared_ptr<IterativeSolver> cost_minimizer_gsl_create(
                  const boost::shared_ptr<CostFunc>& cost, int max_cost_function_calls,
                  double size_tol, const Array<double,1>& init_step_size, bool vrbs)
{
  return boost::shared_ptr<IterativeSolver>(new CostMinimizerGSL( cost, max_cost_function_calls, size_tol, init_step_size, vrbs ));
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(CostMinimizerGSL, IterativeSolver)
.def(luabind::constructor< const boost::shared_ptr<CostFunc>&, int, double, const Array<double,1>&, bool >())
.def(luabind::constructor< const boost::shared_ptr<CostFunc>&, int, double, const Array<double,1>& >())
.def(luabind::constructor< const boost::shared_ptr<CostFunc>&, int, double >())
.def(luabind::constructor< const boost::shared_ptr<CostFunc>&, int >())
.scope
[
 luabind::def("create", &cost_minimizer_gsl_create)
]
REGISTER_LUA_END()
#endif




void print_state( unsigned int iter, gsl_multimin_fminimizer * s, int status)
{
  printf( "Solver '%s';  iter = %3u;  c(x) = %g;  gsl status = %s\n",
          gsl_multimin_fminimizer_name(s), iter,
          gsl_multimin_fminimizer_minimum(s), gsl_strerror((int)status) );

  printf( "Where x is\n" );
  (void) gsl_vector_fprintf(stdout, gsl_multimin_fminimizer_x(s), "%25.14lf");

  printf( "Average distance from the simplex center to its vertices = %25.10lf\n",
          gsl_multimin_fminimizer_size(s) );

  printf( "\n" );
}



CostMinimizerGSL::CostMinimizerGSL(const boost::shared_ptr<CostFunc>& p,
                                   int max_cost_function_calls, double size_tol, 
                                   const Array<double,1>& init_step_size,
                                   bool vrbs)
  : CostMinimizer(p, max_cost_function_calls, vrbs),
    Size_tol(size_tol), Initial_step_size(init_step_size)
{
  if( (init_step_size.size() > 0) &&
      ((int) init_step_size.size() != p->expected_parameter_size())) {
    Exception e;
    e << "If initial-step-size provided, its size must be equal to the expected-parameter-size:\n"
      << " Initial-step-size: " << init_step_size.size() << "\n"
      << " Expected-parameter-size: " << p->expected_parameter_size() << "\n";
    throw e;
  }
}


void CostMinimizerGSL::solve()
{
  // Selecting a problem
  gsl_multimin_function f = gsl_get_mdm(&(*P));

  // select the solver
  const gsl_multimin_fminimizer_type * T = get_gsl_multimin_fminimizer();

  // for gsl status
  int gsl_status = GSL_FAILURE;
  int gsl_status_s = GSL_CONTINUE;

  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, f.n);

  blitz::Array<double, 1> X(P->parameters());

  gsl_vector *ss = gsl_vector_alloc(f.n);
  bool ss_initialized = false;
  //
  if( ss ) {
    //  Choosing a good initial step-size is in itself a
    //  research topic.

    if(Initial_step_size.size() > 0) {
      ss_initialized = (gsl_vector_memcpy(ss, GslVector(Initial_step_size).gsl()) == GSL_SUCCESS);
    } else {

      // Here is a simple choice.
      //
      gsl_vector_set_all(ss, 1.0);

      // Another common choice.
      //
//      for( size_t i=0; i<f.n; i++ )
//        gsl_vector_set( ss, i, (X(i)==0.0)?0.00025:0.05 );

      ss_initialized = true;
    }
  }

  int num_step = 0;
  stat = UNTRIED;
  if( s && ss_initialized )
    if( !(gsl_status = gsl_multimin_fminimizer_set(s, &f, GslVector(X).gsl(), ss)) ) {
      stat = CONTINUE;
      do {
        num_step++;

        // The following two lines are only for recording purpose.
        // They record info at the initial guess (the starting point).
        //
        record_cost_at_accepted_point(P->cost());
        record_accepted_point(P->parameters());

        // A return status of GSL_ENOPROG by the following function
        // does not suggest convergence, and it only means that the
        // function call did not encounter an error.  However, a
        // return status of GSL_ENOPROG means that the function call
        // did not make any progress.  Stalling can happen at a true
        // minimum or some other reason.
        //
        gsl_status = gsl_multimin_fminimizer_iterate(s);
        if( (gsl_status != GSL_SUCCESS) && (gsl_status != GSL_ENOPROG) ) {
          stat = ERROR;
          break;
        }

        gsl_status_s = gsl_multimin_test_size(gsl_multimin_fminimizer_size(s), Size_tol);
        if(gsl_status_s == GSL_SUCCESS) {
          stat = SUCCESS;
          break;
        }

        if(gsl_status == GSL_ENOPROG) {
          stat = STALLED;
          break;
        }

        if( verbose ) print_state(num_step, s, gsl_status_s);
      } while (gsl_status_s == GSL_CONTINUE && P->num_cost_evaluations() < max_cost_f_calls);
    }

  // The following two lines are only for recording purpose.
  //
  record_cost_at_accepted_point(P->cost());
  record_accepted_point(P->parameters());

  if( verbose && (stat != CONTINUE) )
    print_state( num_step, s, ((stat == ERROR)?gsl_status:gsl_status_s) );

  if( ss ) gsl_vector_free(ss);
  if( s ) gsl_multimin_fminimizer_name(s);
}
