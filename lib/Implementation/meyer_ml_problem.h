#ifndef MEYER_ML_PROBLEM
#define MEYER_ML_PROBLEM
#include <max_likelihood.h>
#include <model_measure_meyer.h>


namespace FullPhysics {

class MeyerMLProblem : virtual public MaxLikelihood,
		       virtual public ModelMeasureMeyer {
public:
  MeyerMLProblem(const blitz::Array<double, 1>& measurement, 
                 const blitz::Array<double, 1>& measurement_error_cov);
  virtual ~MeyerMLProblem() {}
  virtual void print(std::ostream& Os) const
  { Os << "MeyerMLProblem"; }
private:
  MeyerMLProblem() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(MeyerMLProblem);
#endif
