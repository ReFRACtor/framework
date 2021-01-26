#ifndef NLLS_MAX_LIKELIHOOD_H
#define NLLS_MAX_LIKELIHOOD_H
#include <max_likelihood.h>
#include <nlls_problem.h>
#include <nlls_problem_state.h>
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/******************************************************************
  This class casts a maximum likelihood problem into the
  Nonlinear Least Squares problem.
*******************************************************************/
class NLLSMaxLikelihood : 
    virtual public NLLSProblem, virtual public NLLSProblemState {

public:

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

  NLLSMaxLikelihood(const boost::shared_ptr<MaxLikelihood>& ml)
    : NLLSProblem(), NLLSProblemState(), ML(ml)
  {}

  virtual ~NLLSMaxLikelihood() {}


//-----------------------------------------------------------------------
/// Return the residual of the NLLS problem
/// at the current set point
/// 
/// \return Residual
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> residual();


//-----------------------------------------------------------------------
/// Return the Jacobian of the residual of the NLLS problem
/// at the current set point.
///
/// \return The Jacobian of the cost function.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> jacobian();


//-----------------------------------------------------------------------
/// Return the size of the residual that will be returned by residual()
//-----------------------------------------------------------------------

  virtual int residual_size() const
  { return ML->measurement_size(); }


//-----------------------------------------------------------------------
/// Return the size of the parameter X.
//-----------------------------------------------------------------------

  virtual int expected_parameter_size() const
  { return ML->expected_parameter_size(); }


//-----------------------------------------------------------------------
/// Sets the problem at a new point in the parameter space.
/// 
/// \param x Input value
//-----------------------------------------------------------------------

  virtual void parameters(const blitz::Array<double, 1>& x);


//-----------------------------------------------------------------------
/// Just returns the current values of parameters.
/// This method is redefined here (see the root base
/// class) because of a compiler bug; otherwise, there
/// should be no need for its redefinition.
/// 
/// \return Current parameter values
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> parameters() const
  { return NLLSProblemState::parameters(); }


  boost::shared_ptr<MaxLikelihood> max_likelihood()
  { return ML; }


//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "NLLSMaxLikelihood"; }


protected:
  boost::shared_ptr<MaxLikelihood> ML;
  NLLSMaxLikelihood() {}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}
FP_EXPORT_KEY(NLLSMaxLikelihood);
#endif
