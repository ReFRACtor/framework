// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "nlls_problem_state.h"
%}
%base_import(problem_state)
%fp_shared_ptr(FullPhysics::NLLSProblemState);

namespace FullPhysics {
class NLLSProblemState : virtual public ProblemState {
public:
  NLLSProblemState();
  NLLSProblemState(const NLLSProblemState& s);
  %pickle_serialization();
};
}
