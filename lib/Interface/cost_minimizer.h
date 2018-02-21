#ifndef COST_MINIMIZER_H
#define COST_MINIMIZER_H
#include <iterative_solver.h>
#include <cost_func.h>


namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for all iterative cost minimizers
///        that do not require derivatives of any order.
///
/// This is the base class for methods that iteratively
/// minimize a scalar cost function.  In other words, this 
/// is a base class for methods that find a point in the
/// parameter space where the cost function is at least 
/// locally minimum.
///
/// This class is associated with a problem (CostFunc)
/// because the problem interface is determined:
///   - provide a point in the parameter space
///   - evaluate the cost function (a scalar function)
///     at the point
//-----------------------------------------------------------------------

class CostMinimizer : 
    public IterativeSolver {

public:


//-----------------------------------------------------------------------
/// \brief Constructor
/// 
/// \param[in] p
///            The cost minimization problem
///
/// \param[in] max_cost_function_calls 
///            read related base class comments
///
/// \param[in] vrbs
///            read related base class comments
//-----------------------------------------------------------------------

  CostMinimizer(const boost::shared_ptr<CostFunc>& p,
                int max_cost_function_calls, bool vrbs)
    : IterativeSolver(max_cost_function_calls, vrbs),
      P(p)
  {}


  virtual ~CostMinimizer() {}


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "CostMinimizer"; }


protected:

  boost::shared_ptr<CostFunc> P;

};
}
#endif
