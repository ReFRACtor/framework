#include "absorber.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Absorber::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableAbsorber);
    
}

FP_IMPLEMENT(Absorber);
FP_OBSERVER_SERIALIZE(Absorber);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Absorber)
.def("gas_index", &Absorber::gas_index)
.def("number_species", &Absorber::number_species)
.def("gas_name", &Absorber::gas_name)
.def("absorber_vmr", &Absorber::absorber_vmr)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Timer for optical_depth_each_layer.
//-----------------------------------------------------------------------

AccumulatedTimer Absorber::timer("Absorber optical_depth_each_layer");

//-----------------------------------------------------------------------
/// Map a gas name to the index number it appears in
/// optical_depth_each_layer. This return -1 if the Name is not one of
/// the gases.
//-----------------------------------------------------------------------

int Absorber::gas_index(const std::string& Name) const
{
  for(int i = 0; i < number_species(); ++i)
    if(gas_name(i) == Name)
      return i;
  return -1;
}


