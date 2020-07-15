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
            temp_levels(Temp_levels), skin_temp_flag(Skin_temp_flag), gas_levels(Gas_levels),
            emissivity_flags(Emissivity_flags), reflectivity_flags(Reflectivity_flags) { }
    virtual ~OssRetrievalFlags() { }
    Array<int, 1> temp_levels; ///< indexes to flag retrieved temperature levels
    bool skin_temp_flag; ///< flag to indicate if skin temperature was retrieved
    std::vector<Array<int, 1>> gas_levels; ///< per gas, indexes to flag retrieved gas levels
    Array<int, 1> emissivity_flags; ///< indexes to flag retrieved emissivities
    Array<int, 1> reflectivity_flags; ///< indexes to flag retrieved reflectivities
    /* TODO: Add support for clouds */
};
}

#endif
