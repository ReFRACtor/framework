#include "pcloud_direct.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(Pcloud,
				 SubStateVectorArrayPcloud);

template<class Archive>
void PcloudDirect::serialize(Archive & ar,
					 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayPcloud)
    & FP_NVP(units);
}

FP_IMPLEMENT(PcloudDirect);
#endif


PcloudDirect::PcloudDirect(const ArrayWithUnit<double, 1>& pcloud, boost::shared_ptr<StateMapping> in_map)
{
    // Save units seperately so it can not be saved into the state vector coefficients
    units = pcloud.units;

    init(pcloud.value, in_map);
}

AutoDerivativeWithUnit<double> PcloudDirect::pressure_cloud(int sensor_index) const
{
  return AutoDerivativeWithUnit<double>(mapped_state()(sensor_index), units);
}

boost::shared_ptr<Pcloud> PcloudDirect::clone() const
{
  return boost::shared_ptr<PcloudDirect>(new PcloudDirect(ArrayWithUnit<double, 1>(coefficient().value(), units),state_mapping()));
}

std::string PcloudDirect::state_vector_name_i(int i) const
{
  std::stringstream sv_name;
  if(mapping->name() != "linear")
    sv_name << mapping->name() << " Pressure Cloud for sensor index " << i;
  else
    sv_name << "Pressure Cloud for sensor index " << i;
  return sv_name.str();
}
