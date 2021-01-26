#ifndef BARD_ML_PROBLEM
#define BARD_ML_PROBLEM
#include <max_likelihood.h>
#include <model_measure_bard.h>


namespace FullPhysics {

class BardMLProblem : virtual public MaxLikelihood,
		      virtual public ModelMeasureBard {
public:
  BardMLProblem(const blitz::Array<double, 1>& measurement, 
                const blitz::Array<double, 1>& measurement_error_cov);
  virtual ~BardMLProblem() {}
  virtual void print(std::ostream& Os) const
  { Os << "BardMLProblem"; }
private:
  BardMLProblem() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(BardMLProblem);
#endif
