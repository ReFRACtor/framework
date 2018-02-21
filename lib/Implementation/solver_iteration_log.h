#ifndef SOLVER_ITER_LOG_H
#define SOLVER_ITER_LOG_H
#include "iterative_solver.h"
#include "state_vector.h"
#include "logger.h"

namespace FullPhysics {
  /****************************************************************//**
    This is an observer of IterativeSolver that writes out the 
    state vector values in a nicely formated way on each accepted
    iteration.

    For this to work it must be registered as an Observer of the
    instantiated IterativeSolver.

  *******************************************************************/

  class SolverIterationLog : public Printable<SolverIterationLog>, public Observer<IterativeSolver> {
  public:
      static const int SV_PRINT_WIDTH = 17;
    
      SolverIterationLog(const boost::shared_ptr<StateVector>& Sv) : sv_obj(Sv) {};

      void notify_update(const IterativeSolver& solver);

      virtual void print(std::ostream& Os) const { Os << "SolverIterationLog";}
  private:
      boost::shared_ptr<StateVector> sv_obj;
};
}
#endif
