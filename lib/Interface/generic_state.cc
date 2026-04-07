#include "generic_state.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void GenericState::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableGenericState);
    
}

SUB_STATE_VECTOR_ARRAY_SERIALIZE(GenericState,
				 SubStateVectorArrayGenericState);

template<class Archive>
void GenericStateImpBase::serialize(Archive & ar,
					 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorArrayGenericState);
}

FP_IMPLEMENT(GenericStateImpBase);
FP_IMPLEMENT(GenericState);
FP_OBSERVER_SERIALIZE(GenericState);
#endif

std::string GenericStateImpBase::state_vector_name_i(int i) const
{
  std::stringstream sv_name;
  if(mapping->name() != "linear")
    sv_name << mapping->name() << " Generic State " << i;
  else
    sv_name << "Generic State " << i;
  return sv_name.str();
}
