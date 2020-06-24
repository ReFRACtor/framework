#include <model_measure_bard.h>
#include "fp_serialize_support.h"


using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ModelMeasureBard::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ModelMeasure);
}

FP_IMPLEMENT(ModelMeasureBard);
#endif


void ModelMeasureBard::model_eval()
{
  if(M.size() <= 0) {
    assert_parameter_set_correctly();
    M.resize(measurement_size());

    for(int i=1; i<=measurement_size(); i++)
      M(i-1) = X(0) + i/((16-i)*X(1) + (8-abs(8-i))*X(2));
  }
}



void ModelMeasureBard::jacobian_eval()
{
  if(K.size() <= 0) {
    assert_parameter_set_correctly();
    K.resize(measurement_size(), expected_parameter_size());

    for(int i=1; i<=measurement_size(); i++) {
      double denom = (16-i)*X(1) + (8-abs(8-i))*X(2);
      denom = denom*denom;
      K(i-1,0) = 1.0;  K(i-1,1) = -(16-i)*i/denom;  K(i-1,2) = -(8-abs(8-i))*i/denom;
    }
  }
}


void ModelMeasureBard::model_jacobian_eval()
{
  model_eval();
  jacobian_eval();
}
