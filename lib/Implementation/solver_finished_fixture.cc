#include "solver_finished_fixture.h"
#include "ifstream_cs.h"

using namespace FullPhysics;

SolverFinishedFixture::SolverFinishedFixture()
{
  solver = config_solver;
  IfstreamCs in(test_data_dir() + "in/solver/connor_converged");
  in >> *solver;
  config_state_vector->update_state(solver->x_solution(), 
				    solver->aposteriori_covariance());
  initial_sv.reference(config_initial_guess->initial_guess());
  apriori_sv.reference(config_initial_guess->apriori());
  apriori_cov.reference(config_initial_guess->apriori_covariance());

  error_analysis.reset(new ErrorAnalysis(solver, config_atmosphere, config_forward_model, config_observation));
}
