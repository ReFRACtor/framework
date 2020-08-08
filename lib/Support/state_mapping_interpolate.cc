#include "state_mapping_interpolate.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<template<class TX, class TY> class Interp> template<class Archive>
void StateMappingInterpolate<Interp>::serialize(Archive& ar, const unsigned int UNUSED(version))
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateMapping)
    & FP_NVP(press_to) & FP_NVP(press_from) & FP_NVP(map_name);
}

//FP_IMPLEMENT(StateMappingInterpolateLinearLinear);
//FP_IMPLEMENT(StateMappingInterpolateLogLinear);
//FP_IMPLEMENT(StateMappingInterpolateLogLog);
//FP_IMPLEMENT(StateMappingInterpolateLinearLog);
#endif
