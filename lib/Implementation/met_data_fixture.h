#ifndef MET_DATA_FIXTURE
#define MET_DATA_FIXTURE

#include "configuration_fixture.h"
#include "unit_test_support.h"
#include "meteorology.h"

using namespace FullPhysics;
using namespace blitz;

class MetDataFixture : public GlobalFixture {
public:
    boost::shared_ptr<Meteorology> met_data;
    blitz::Array<double, 1> pressure_in;
    double psurf_in;
    MetDataFixture();
};

#endif
