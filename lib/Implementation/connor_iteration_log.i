%include "fp_common.i"
%{
#include "connor_iteration_log.h"
%}

%base_import(connor_solver)
%import "state_vector.i"

%fp_shared_ptr(FullPhysics::ConnorIterationLog);

%template(ObserverConnorIterationLog) FullPhysics::Observer<FullPhysics::ConnorSolver>;

namespace FullPhysics {
  class ConnorIterationLog : public Observer<ConnorSolver> {
  public:
    ConnorIterationLog(const boost::shared_ptr<StateVector>& Sv);
    void notify_update(const ConnorSolver& solver);
    std::string print_to_string() const;
};
}
