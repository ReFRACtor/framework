#include <bard_ml_problem.h>
#include "fp_serialize_support.h"


using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void BardMLProblem::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ModelMeasureBard)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MaxLikelihood);
}

FP_IMPLEMENT(BardMLProblem);
#endif

BardMLProblem::BardMLProblem(
    const blitz::Array<double, 1>& measurement, 
    const blitz::Array<double, 1>& measurement_error_cov)
  : ModelMeasure(measurement, measurement_error_cov),
    ModelMeasureBard()
{
}
