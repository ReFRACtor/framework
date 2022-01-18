#ifndef SERIALIZED_CONFIGURATION_FIXTURE_H
#define SERIALIZED_CONFIGURATION_FIXTURE_H

#include "configuration_fixture.h"
#include "fp_exception.h"

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

}
#endif
