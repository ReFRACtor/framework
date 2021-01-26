#include "surface_temperature_direct.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

SUB_STATE_VECTOR_ARRAY_SERIALIZE(SurfaceTemperature,
				 SubStateVectorArraySurfaceTemperature);

template<class Archive>
void SurfaceTemperatureDirect::serialize(Archive & ar,
					 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArraySurfaceTemperature)
    & FP_NVP(units);
}

FP_IMPLEMENT(SurfaceTemperatureDirect);
#endif


SurfaceTemperatureDirect::SurfaceTemperatureDirect(const ArrayWithUnit<double, 1>& surf_temp)
{
    // Save units seperately so it can not be saved into the state vector coefficients
    units = surf_temp.units;

    init(surf_temp.value);
}

AutoDerivativeWithUnit<double> SurfaceTemperatureDirect::surface_temperature(int channel_index) const
{
    return AutoDerivativeWithUnit<double>(coeff(channel_index), units);
}

boost::shared_ptr<SurfaceTemperature> SurfaceTemperatureDirect::clone() const
{
    return boost::shared_ptr<SurfaceTemperatureDirect>(new SurfaceTemperatureDirect(ArrayWithUnit<double, 1>(coefficient().value(), units)));
}

std::string SurfaceTemperatureDirect::state_vector_name_i(int i) const
{
    std::stringstream sv_name;
    sv_name << "Surface Temperature for channel " << i;
    return sv_name.str();
}
