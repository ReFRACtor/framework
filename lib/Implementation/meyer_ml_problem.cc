#include <meyer_ml_problem.h>
#include "fp_serialize_support.h"


using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void MeyerMLProblem::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ModelMeasureMeyer)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MaxLikelihood);
}

FP_IMPLEMENT(MeyerMLProblem);
#endif

MeyerMLProblem::MeyerMLProblem(
    const blitz::Array<double, 1>& measurement, 
    const blitz::Array<double, 1>& measurement_error_cov)
  : ModelMeasure(measurement, measurement_error_cov),
    ModelMeasureMeyer()
{
}
