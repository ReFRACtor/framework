#include "sample_grid.h"
using namespace FullPhysics;
using namespace blitz;

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

