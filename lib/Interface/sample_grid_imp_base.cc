#include "sample_grid_imp_base.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(SampleGrid, SubStateVectorArraySampleGrid);

template<class Archive>
void SampleGridImpBase::serialize(Archive & ar,
                                  const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArraySampleGrid);
}

FP_IMPLEMENT(SampleGridImpBase);
#endif
