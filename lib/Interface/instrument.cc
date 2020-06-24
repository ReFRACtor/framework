#include "instrument.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Instrument::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableInstrument);
    
}

FP_IMPLEMENT(Instrument);
FP_OBSERVER_SERIALIZE(Instrument);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Instrument)
.def("number_spectrometer", &Instrument::number_spectrometer)
REGISTER_LUA_END()
#endif
