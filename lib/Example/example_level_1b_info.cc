#include "example_level_1b_info.h"

using namespace FullPhysics;
using namespace blitz;

ExampleLevel1bInfo::ExampleLevel1bInfo(const boost::shared_ptr<HdfFile>& input_file)
:input(input_file)
{
    // All work done in init
}

ExampleLevel1bInfo::ExampleLevel1bInfo(const std::string& input_filename)
:input(new HdfFile(input_filename))
{
    // All work done in init    
}

std::vector<std::shared_ptr<Level1b>> ExampleLevel1bInfo::level1b_list() {
    std::string dataset_name = "/observation_ids";
    TinyVector<int, 1> ds_shape = input->read_shape<1>(dataset_name);
    TinyVector<int, 1> start, size;
    start = 0;
    size = ds_shape(0);
    Array<std::string, 1> obs_arr = input->read_field<std::string, 1>(dataset_name, start, size);

    std::vector<std::shared_ptr<Level1b>> level1b_series = std::vector<std::shared_ptr<Level1b>>();
    for (auto obs_str : obs_arr) {
        std::shared_ptr<Level1b> new_level1b_ptr (new ExampleLevel1b(input, obs_str));
        level1b_series.push_back(new_level1b_ptr);
    }

    return level1b_series;

}
