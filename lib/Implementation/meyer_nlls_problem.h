#ifndef MEYER_NLLS_PROBLEM
#define MEYER_NLLS_PROBLEM
#include "nlls_problem.h"
#include "nlls_problem_state.h"


namespace FullPhysics {

class MeyerNLLSProblem : virtual public NLLSProblem,
			 virtual public NLLSProblemState {
public:
  MeyerNLLSProblem()
    : NLLSProblem()
  {}
  virtual ~MeyerNLLSProblem() {}
  virtual int residual_size() const { return 16; }
  virtual int expected_parameter_size() const { return 3; }
  virtual blitz::Array<double, 1> residual();
  virtual blitz::Array<double, 2> jacobian();
  virtual void print(std::ostream& Os) const
    { Os << "MeyerNLLSProblem"; }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(MeyerNLLSProblem);
#endif
