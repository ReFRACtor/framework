#ifndef BARD_NLLS_PROBLEM
#define BARD_NLLS_PROBLEM
#include "nlls_problem.h"
#include "nlls_problem_state.h"


namespace FullPhysics {

class BardNLLSProblem : virtual public NLLSProblem,
			virtual public NLLSProblemState {
public:
  BardNLLSProblem()
    : NLLSProblem()
  {}
  virtual ~BardNLLSProblem() {}
  virtual int residual_size() const { return 15; }
  virtual int expected_parameter_size() const { return 3; }
  virtual blitz::Array<double, 1> residual();
  virtual blitz::Array<double, 2> jacobian();
  virtual void print(std::ostream& Os) const
    { Os << "BardNLLSProblem"; }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(BardNLLSProblem);
#endif
