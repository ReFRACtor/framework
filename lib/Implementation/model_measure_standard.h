#ifndef MODEL_MEASURE_STANDARD_H
#define MODEL_MEASURE_STANDARD_H
#include "model_measure.h"
#include "forward_model.h"
#include "observation.h"
#include "state_vector.h"

#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/******************************************************************
  This class implements what is common to 
    - "maximum a posteriori" for a standard forward model, and
    - "maximum likelihood" for a standard forward model
*******************************************************************/
class ModelMeasureStandard : 
  virtual public ModelMeasure {

public:

  virtual ~ModelMeasureStandard() {}

  virtual void model_eval();

  virtual void jacobian_eval();

  virtual void model_jacobian_eval();

  virtual int expected_parameter_size() const;


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
  { return ModelMeasure::parameters(); }


//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "ModelMeasureStandard"; }

//-----------------------------------------------------------------------
/// Underlying forward model.
//-----------------------------------------------------------------------

  const std::vector<boost::shared_ptr<ForwardModel> >& forward_model() const
  { return fm; }

//-----------------------------------------------------------------------
/// Underlying observation.
//-----------------------------------------------------------------------

  const std::vector<boost::shared_ptr<Observation> >& observation() const
  { return obs; }

//-----------------------------------------------------------------------
/// Underlying 
//-----------------------------------------------------------------------

  const boost::shared_ptr<StateVector>& state_vector() const { return sv; };

  virtual void measurement_eval();
  virtual void measurement_jacobian_eval();
protected:

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

  ModelMeasureStandard(const boost::shared_ptr<ForwardModel>& forward_model, const boost::shared_ptr<Observation>& observation, const boost::shared_ptr<StateVector>& state_vector);

  ModelMeasureStandard(const std::vector<boost::shared_ptr<ForwardModel> >& forward_model, const std::vector<boost::shared_ptr<Observation> >& observation, const boost::shared_ptr<StateVector>& state_vector);
  
  ModelMeasureStandard() {}

  virtual void radiance_from_fm(bool skip_check=false);

  std::vector<boost::shared_ptr<ForwardModel> > fm;
  std::vector<boost::shared_ptr<Observation> > obs;
  boost::shared_ptr<StateVector> sv;
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(ModelMeasureStandard);
#endif
