%include "fp_common.i"
%{
#include "solver_iteration_log.h"
%}

%base_import(iterative_solver)
%import "state_vector.i"

%fp_shared_ptr(FullPhysics::SolverIterationLog);

%template(ObserverSolverIterationLog) FullPhysics::Observer<FullPhysics::IterativeSolver>;

namespace FullPhysics {
  class SolverIterationLog : public Observer<IterativeSolver> {
  public:
    SolverIterationLog(const boost::shared_ptr<StateVector>& Sv);
    void notify_update(const IterativeSolver& solver);
    std::string print_to_string() const;
};
}
