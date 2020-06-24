%include "fp_common.i"
%{
#include "iterative_solver.h"
%}

%base_import(observer)

%import "observer.i"

%fp_shared_ptr(FullPhysics::IterativeSolver);

%fp_shared_ptr(FullPhysics::Observable<FullPhysics::IterativeSolver>);
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::IterativeSolver>);
%template(ObservableIterativeSolver) FullPhysics::Observable<FullPhysics::IterativeSolver>;

// It is very important that the dirctor come before the template or else the director won't work properly 
%feature("director") FullPhysics::Observer<FullPhysics::IterativeSolver>;
%template(ObserverIterativeSolver) FullPhysics::Observer<FullPhysics::IterativeSolver>;

namespace FullPhysics {
class IterativeSolver : public Observable<IterativeSolver> {
public:
  enum status_t {SUCCESS, CONTINUE, STALLED, ERROR, UNTRIED};
  IterativeSolver(int max_cost_function_calls, bool vrbs);
  virtual ~IterativeSolver();
  virtual void add_observer(Observer<IterativeSolver>& Obs);
  virtual void remove_observer(Observer<IterativeSolver>& Obs);
  %python_attribute(num_accepted_steps, int)
  %python_attribute(accepted_points, std::vector< blitz::Array<double, 1> >)
  %python_attribute(cost_at_accepted_points, std::vector<double>)
  virtual void solve() = 0;
  %python_attribute(status, status_t)
  %python_attribute(status_str, std::string)
  std::string print_to_string() const;
  %pickle_serialization();
};
}
