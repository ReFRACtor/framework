#include "scale_cloud_direct.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(ScaleCloud,
				 SubStateVectorArrayScaleCloud);

template<class Archive>
void ScaleCloudDirect::serialize(Archive & ar,
					 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayScaleCloud);
}

FP_IMPLEMENT(ScaleCloudDirect);
#endif


ScaleCloudDirect::ScaleCloudDirect(const blitz::Array<double, 1>& scale_cloud, boost::shared_ptr<StateMapping> in_map)
{
  init(scale_cloud, in_map);
}

AutoDerivative<double> ScaleCloudDirect::scale_cloud(int sensor_index) const
{
  return mapped_state()(sensor_index);
}

boost::shared_ptr<ScaleCloud> ScaleCloudDirect::clone() const
{
  return boost::shared_ptr<ScaleCloudDirect>(new ScaleCloudDirect(coefficient().value(), state_mapping()));
}

std::string ScaleCloudDirect::state_vector_name_i(int i) const
{
  std::stringstream sv_name;
  if(mapping->name() != "linear")
    sv_name << mapping->name() << " Scale Cloud for sensor index " << i;
  else
    sv_name << "Scale Cloud for sensor index " << i;
  return sv_name.str();
}
