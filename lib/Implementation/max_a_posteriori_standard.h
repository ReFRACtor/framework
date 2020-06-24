#ifndef MAX_A_POSTERIORI_STANDARD_H
#define MAX_A_POSTERIORI_STANDARD_H
#include <max_a_posteriori.h>
#include <model_measure_standard.h>

namespace FullPhysics {
/******************************************************************
  This class implements "maximum a posteriori" for a standard forward model
*******************************************************************/
class MaxAPosterioriStandard : 
  virtual public MaxAPosteriori, virtual public ModelMeasureStandard {

public:

  MaxAPosterioriStandard(const boost::shared_ptr<ForwardModel>& fm,
          const boost::shared_ptr<Observation>& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> a_priori_cov);

  MaxAPosterioriStandard(const std::vector<boost::shared_ptr<ForwardModel> >& fm,
 	  const std::vector<boost::shared_ptr<Observation> >& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> a_priori_cov);
  
  virtual ~MaxAPosterioriStandard() {}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "MaxAPosterioriStandard"; }

protected:

  //  TEMPORARY
  //
  // Should go away after we end support for 
  // fixed pressure level grid.
  virtual void vanishing_params_update();
private:
  MaxAPosterioriStandard() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(MaxAPosterioriStandard);
#endif
