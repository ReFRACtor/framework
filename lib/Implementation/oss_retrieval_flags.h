#ifndef OSS_RETRIEVAL_FLAGS_H
#define OSS_RETRIEVAL_FLAGS_H

#include <vector>
#include <blitz/array.h>

#include "generic_object.h"

namespace FullPhysics {
/****************************************************************//**
  This class tracks the OSS FM retrieval flags / retrieved levels
*******************************************************************/

class OssRetrievalFlags : public GenericObject {
public:
    OssRetrievalFlags(const blitz::Array<int, 1>& Temp_levels, 
                      const blitz::Array<bool, 1>& Skin_temp_sensors,
                      const std::vector<blitz::Array<int, 1> >& Gas_levels,
                      const blitz::Array<int, 1>& Emissivity_flags,
                      const blitz::Array<int, 1>& Reflectivity_flags) 

        : temp_levels_(Temp_levels), skin_temp_sensors_(Skin_temp_sensors), gas_levels_(Gas_levels),
          emissivity_flags_(Emissivity_flags), reflectivity_flags_(Reflectivity_flags) { }

    virtual ~OssRetrievalFlags() { }

    const blitz::Array<int, 1> temp_levels() const { return temp_levels_; }
    const blitz::Array<bool, 1> skin_temp_sensors() const { return skin_temp_sensors_; }
    const std::vector<blitz::Array<int, 1>> gas_levels() const { return gas_levels_; }
    const blitz::Array<int, 1> emissivity_flags() const { return emissivity_flags_; }
    const blitz::Array<int, 1> reflectivity_flags() const { return reflectivity_flags_; }

    int num_total_flags() const;

private:

    blitz::Array<int, 1> temp_levels_; ///< indexes to flag retrieved temperature levels
    blitz::Array<bool, 1> skin_temp_sensors_; ///< flag to indicate if skin temperature for each sensor (band)
    std::vector<blitz::Array<int, 1>> gas_levels_; ///< per gas, indexes to flag retrieved gas levels
    blitz::Array<int, 1> emissivity_flags_; ///< indexes to flag retrieved emissivities
    blitz::Array<int, 1> reflectivity_flags_; ///< indexes to flag retrieved reflectivities

    /* TODO: Add support for clouds */
};
}

#endif
