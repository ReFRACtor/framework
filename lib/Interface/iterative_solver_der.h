#ifndef ITERATIVE_SOLVER_DER_H
#define ITERATIVE_SOLVER_DER_H
#include <iterative_solver.h>
#include <cost_func_diff.h>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for all iterative optimizers that use
///        first order derivatives.
///
/// This class is the base class for iterative optimizers that 
/// use first order derivatives.
///
/// Similar to its base class IterativeSolver, IterativeSolverDer
/// is also not associated with any problem for the same reason
/// mentioned in the comment section of IterativeSolver class.
//-----------------------------------------------------------------------

class IterativeSolverDer : 
    public IterativeSolver {

public:


//-----------------------------------------------------------------------
/// \brief Constructor
/// 
/// \param[in] max_cost_function_calls 
///            read related base class comments
///
/// \param[in] vrbs
///            read related base class comments
//-----------------------------------------------------------------------

  IterativeSolverDer(int max_cost_function_calls, bool vrbs)
    : IterativeSolver(max_cost_function_calls, vrbs)
  {}


  virtual ~IterativeSolverDer() {}


//-----------------------------------------------------------------------
/// \brief Returns a vector (std) of gradients evaluated
///        at accepted points.
///
/// This method returns a std vector of gradients
/// computed at the accepted points.  In other words, if the 
/// accepted points and the computed gradients at these points
/// are recorded correctly, then
///   - gradient_at_accepted_points()[0] is the gradient of the cost 
///     function evaluated at accepted_points()[0]
///   - gradient_at_accepted_points()[1] is the gradient of the cost
///     function evaluated at accepted_points()[1]
///   - ...
///   - and finally gradient_at_accepted_points()[num_accepted_steps()]
///     is the gradient of the cost function evaluated at 
///     accepted_points()[num_accepted_steps()]
//-----------------------------------------------------------------------

  virtual std::vector< blitz::Array<double, 1> > gradient_at_accepted_points() const
  { return Gradient_at_accepted_points; }


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "IterativeSolverDer"; }


protected:

//-----------------------------------------------------------------------
/// \brief For recording the gradient of the cost 
///        function evaluated at an accepted point
///
/// This method is called to record the gradient of the cost function
/// evaluated at an accepted point.  It is the responsibility of the
/// implementer of the solve() method to record the gradients evaluated
/// at the accepted points.  The gradients must be recorded in the same
/// order that they are evaluated.
///
/// \param[in] gradient
///            gradient of the cost function evaluated at
///            an accepted point in the parameter space
//-----------------------------------------------------------------------

  void record_gradient_at_accepted_point(const blitz::Array<double, 1>& gradient)
  { Gradient_at_accepted_points.push_back(gradient); }


  IterativeSolverDer() {}
private:

  std::vector< blitz::Array<double, 1> > Gradient_at_accepted_points;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(IterativeSolverDer);
#endif
