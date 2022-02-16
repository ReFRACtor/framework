#include "fluorescence_effect.h"
#include "serialized_configuration_fixture.h"
#include <boost/shared_ptr.hpp>

using namespace FullPhysics;
using namespace blitz;

class FluorescenceFixture: public SerializedConfigurationFixture {
public:
    FluorescenceFixture() : SerializedConfigurationFixture("fluorescence_example_config.bin.gz")
    {
        config_fluor = obj_map->get<FluorescenceEffect>("fluorescence");
      
        if (!config_fluor) {
            throw Exception("Could not cast SpectrumEffect into FluorescenceEffect");
        }
    }

    virtual ~FluorescenceFixture() {}
    boost::shared_ptr<FluorescenceEffect> config_fluor;
};
