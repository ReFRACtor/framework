// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "solver_iteration_log.h"
%}

%base_import(observer)
%import "iterative_solver.i"
%import "state_vector.i"

%fp_shared_ptr(FullPhysics::SolverIterationLog);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::IterativeSolver>);

%template(ObserverSolverIterationLog) FullPhysics::Observer<FullPhysics::IterativeSolver>;

namespace FullPhysics {
  class SolverIterationLog : public Observer<IterativeSolver> {
  public:
    SolverIterationLog(const boost::shared_ptr<StateVector>& Sv);
    void notify_update(const IterativeSolver& solver);
    std::string print_to_string() const;
};
}
