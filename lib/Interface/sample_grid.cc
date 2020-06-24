#include "sample_grid.h"
#include "fp_serialize_support.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SampleGrid::serialize(Archive& ar,
			 const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateVectorObserver)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ObservableSampleGrid);
    
}

FP_IMPLEMENT(SampleGrid);
FP_OBSERVER_SERIALIZE(SampleGrid);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_CLASS(SampleGrid)
.def("sample_grid", &SampleGrid::sample_grid)
.def("pixel_grid", &SampleGrid::pixel_grid)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototypes
typedef void(std::vector<boost::shared_ptr<SampleGrid> >::*pbt1)(
        const std::vector<boost::shared_ptr<SampleGrid> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<SampleGrid> >, VectorSampleGrid)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<SampleGrid> >::push_back))
REGISTER_LUA_END()

#endif

