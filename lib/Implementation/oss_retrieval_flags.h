#ifndef OSS_RETRIEVAL_FLAGS_H
#define OSS_RETRIEVAL_FLAGS_H

#include <vector>
#include <blitz/array.h>

#include "generic_object.h"

using namespace FullPhysics;
using namespace blitz;

namespace FullPhysics {
/****************************************************************//**
  This class tracks the OSS FM retrieval flags / retrieved levels
*******************************************************************/

class OssRetrievalFlags : public GenericObject {
public:
    OssRetrievalFlags(Array<int, 1> Temp_levels, bool Skin_temp_flag,
            std::vector<Array<int, 1>> Gas_levels, Array<int, 1> Emissivity_flags,
            Array<int, 1> Reflectivity_flags) :
            temp_levels_(Temp_levels), skin_temp_flag_(Skin_temp_flag), gas_levels_(Gas_levels),
            emissivity_flags_(Emissivity_flags), reflectivity_flags_(Reflectivity_flags) { }
    virtual ~OssRetrievalFlags() { }

    const Array<int, 1> temp_levels() const { return temp_levels_; }
    const bool skin_temp_flag() const { return skin_temp_flag_; }
    const std::vector<Array<int, 1>> gas_levels() const { return gas_levels_; }
    const Array<int, 1> emissivity_flags() const { return emissivity_flags_; }
    const Array<int, 1> reflectivity_flags() const { return reflectivity_flags_; }

    int num_total_flags() const;

private:
    Array<int, 1> temp_levels_; ///< indexes to flag retrieved temperature levels
    bool skin_temp_flag_; ///< flag to indicate if skin temperature was retrieved
    std::vector<Array<int, 1>> gas_levels_; ///< per gas, indexes to flag retrieved gas levels
    Array<int, 1> emissivity_flags_; ///< indexes to flag retrieved emissivities
    Array<int, 1> reflectivity_flags_; ///< indexes to flag retrieved reflectivities

    /* TODO: Add support for clouds */
};
}

#endif
