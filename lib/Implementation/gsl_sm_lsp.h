#ifndef GSL_SM_LSP_H
#define GSL_SM_LSP_H


#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <nlls_problem.h>


int gsl_sm_lsp_f(const gsl_vector *x, void *data, gsl_vector *f);
int gsl_sm_lsp_j(const gsl_vector *x, void *data, gsl_matrix *j);
gsl_multifit_nlinear_fdf gsl_sm_get_lsp_fdf(const FullPhysics::NLLSProblem *lsp);


#endif  /* GSL_SM_LSP_H */
