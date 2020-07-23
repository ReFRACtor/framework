#include "oss_retrieval_flags.h"

using namespace FullPhysics;
using namespace blitz;

int OssRetrievalFlags::num_total_flags() const
{
    int num_flags =  0;

    num_flags += temp_levels_.rows();

    num_flags += blitz::sum(skin_temp_sensors_);

    for (const Array<int, 1>& gas_level : gas_levels_) {
        if(gas_level.rows()) {
            num_flags += gas_level.rows();
        }
    }

    num_flags += emissivity_flags_.rows();
    num_flags += reflectivity_flags_.rows();

    return num_flags;
}
