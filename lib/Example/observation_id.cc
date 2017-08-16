#include "observation_id.h"

using namespace FullPhysics;
using namespace blitz;

template<class T>
ObservationId<T>::ObservationId(const boost::shared_ptr<HdfFile>& input, const std::string& dataset_name, const T& identifier)
{
    Array<T, 1> file_ids = input->read_field<T, 1>(dataset_name);

    data_index_ = -1;
    bool found = false;
    for (int idx = 0; !found && idx < file_ids.rows(); idx++) {
        if(file_ids(idx) == identifier) {
            data_index_ = idx;
            found = true;
        }
    }

    if(!found) {
        Exception err;
        err << "The identifier " << identifier << " was not found in dataset: " << dataset_name << " in file: " << input->file_name();
        throw err;
    }
}

template class ObservationId<std::string>;
template class ObservationId<int>;
