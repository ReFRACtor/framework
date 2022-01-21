#ifndef SERIALIZED_CONFIGURATION_FIXTURE_H
#define SERIALIZED_CONFIGURATION_FIXTURE_H

#include "configuration_fixture.h"
#include "fp_exception.h"
#include "generic_object_map.h"

namespace FullPhysics {

/****************************************************************//**
 This fixture loads configuration objects from a Boost Serialization
 file. This file should have used a GenericObjectMap to store the
 objects named in ConfigurationFixture
*******************************************************************/
class SerializedConfigurationFixture: public ConfigurationFixture {
public:

    SerializedConfigurationFixture(const std::string& Serialized_file);

    virtual ~SerializedConfigurationFixture() = default;

protected:

    virtual void init_variables();
    virtual void init_epsilon();

    boost::shared_ptr<GenericObjectMap> obj_map;

private:

    std::string serialized_filename;
    blitz::Array<double, 1> sv_initial;

};

class LambertianConfigurationFixture : public SerializedConfigurationFixture {
public:
    LambertianConfigurationFixture() : SerializedConfigurationFixture("lambertian_example_config.bin.gz") { ; }
};

class CoxmunkConfigurationFixture : public SerializedConfigurationFixture {
public:
    CoxmunkConfigurationFixture() : SerializedConfigurationFixture("coxmunk_example_config.bin.gz") { ; }
};

class CoxmunkPlusLambertianConfigurationFixture : public SerializedConfigurationFixture {
public:
    CoxmunkPlusLambertianConfigurationFixture() : SerializedConfigurationFixture("coxmunk_lambertian_example_config.bin.gz") { ; }
};

class BrdfVegConfigurationFixture : public SerializedConfigurationFixture {
public:
    BrdfVegConfigurationFixture() : SerializedConfigurationFixture("brdf_veg_example_config.bin.gz") { ; }
};

class BrdfSoilConfigurationFixture : public SerializedConfigurationFixture {
public:
    BrdfSoilConfigurationFixture() : SerializedConfigurationFixture("brdf_soil_example_config.bin.gz") { ; }
};

class TwoBroadenerConfigurationFixture : public SerializedConfigurationFixture {
public:
    TwoBroadenerConfigurationFixture() : SerializedConfigurationFixture("two_broadener_example_config.bin.gz") { ; }
};

}
#endif
