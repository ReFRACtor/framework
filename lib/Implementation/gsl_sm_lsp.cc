#include <stdlib.h>
#include <fp_gsl_matrix.h>
#include <gsl_sm_lsp.h>


using namespace FullPhysics;


int gsl_sm_lsp_f(const gsl_vector *x, void *data, gsl_vector *f)
{
  FullPhysics::NLLSProblem *lsp = (FullPhysics::NLLSProblem *) data;
  blitz::Array<double, 1> b_x(GslVector(const_cast<gsl_vector*>(x), false).blitz_array());
  blitz::Array<double, 1> b_f(lsp->residual_x(b_x));
  gsl_vector_memcpy(f, GslVector(b_f).gsl());

  return GSL_SUCCESS;
}

int gsl_sm_lsp_j(const gsl_vector *x, void *data, gsl_matrix *j)
{
  FullPhysics::NLLSProblem *lsp = (FullPhysics::NLLSProblem *) data;
  blitz::Array<double, 1> b_x(GslVector(const_cast<gsl_vector*>(x), false).blitz_array());
  blitz::Array<double, 2> b_j(lsp->jacobian_x(b_x));
  gsl_matrix_memcpy(j, GslMatrix(b_j).gsl());

  return GSL_SUCCESS;
}

gsl_multifit_nlinear_fdf gsl_sm_get_lsp_fdf(const FullPhysics::NLLSProblem *lsp)
{
  gsl_multifit_nlinear_fdf fdf;
  fdf.f = gsl_sm_lsp_f;
  fdf.df = gsl_sm_lsp_j;  /* set to NULL for finite-difference Jacobian */
  fdf.fvv = NULL;         /* no geodesic acceleration */
  fdf.n = lsp->residual_size();
  fdf.p = lsp->expected_parameter_size();
  fdf.params = (void *) lsp;
  return fdf;
}
