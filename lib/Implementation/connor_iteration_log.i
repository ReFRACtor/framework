%include "common.i"
%{
#include "connor_iteration_log.h"
%}

%base_import(observer)
%import "connor_solver.i"
%import "state_vector.i"

%fp_shared_ptr(FullPhysics::ConnorIterationLog);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::ConnorSolver>);

%template(ObserverConnorIterationLog) FullPhysics::Observer<FullPhysics::ConnorSolver>;

namespace FullPhysics {
  class ConnorIterationLog : public Observer<ConnorSolver> {
  public:
    ConnorIterationLog(const boost::shared_ptr<StateVector>& Sv);
    void notify_update(const ConnorSolver& solver);
    std::string print_to_string() const;
};
}
