#ifndef MET_DATA_FIXTURE
#define MET_DATA_FIXTURE

#include "global_fixture.h"
#include "unit_test_support.h"
#include "meteorology.h"
#include "pressure.h"

using namespace FullPhysics;
using namespace blitz;

class MetDataFixture : public GlobalFixture {
public:
    boost::shared_ptr<Meteorology> met_data;
    boost::shared_ptr<Pressure> pressure;
    boost::shared_ptr<Pressure> pressure_reverse;
    MetDataFixture();
};

#endif
