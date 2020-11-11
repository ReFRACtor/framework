#ifndef COST_MINIMIZER_GSL_H
#define COST_MINIMIZER_GSL_H
#include <gsl/gsl_multimin.h>
#include <cost_minimizer.h>


namespace FullPhysics {
/******************************************************************
  This class is the base class for cost function minimizers
  based on the GSL library.
*******************************************************************/
class CostMinimizerGSL : 
    public CostMinimizer {

public:

//-----------------------------------------------------------------------
/// Initializes the minimizer.
/// 
/// \param p Input value
/// \param max_cost_function_calls Input value
/// \param size_tol
/// \param init_step_size The initial step stize
/// \param vrbs Input value
//-----------------------------------------------------------------------

  CostMinimizerGSL(const boost::shared_ptr<CostFunc>& p,
                   int max_cost_function_calls, double size_tol=0.001, 
                   const blitz::Array<double,1>& init_step_size=blitz::Array<double,1>(),
                   bool vrbs=false);

  virtual ~CostMinimizerGSL() {}

  virtual void solve();

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "CostMinimizerGSL"; }

protected:

  virtual const gsl_multimin_fminimizer_type* get_gsl_multimin_fminimizer()
  { return gsl_multimin_fminimizer_nmsimplex2; /*default*/ }

  double Size_tol;
  blitz::Array<double, 1> Initial_step_size;
  CostMinimizerGSL() {}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(CostMinimizerGSL);
#endif
