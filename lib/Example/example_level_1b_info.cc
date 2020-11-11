#include "example_level_1b_info.h"
#include "fp_serialize_support.h"
#include <boost/make_shared.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void ExampleLevel1bInfo::serialize(Archive & ar,
			       const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Level1bInfo)
    & FP_NVP_(input);
}

FP_IMPLEMENT(ExampleLevel1bInfo);
#endif

ExampleLevel1bInfo::ExampleLevel1bInfo
(const boost::shared_ptr<HdfFile>& input_file)
  :input_(input_file)
{
    // All work done in init
}

ExampleLevel1bInfo::ExampleLevel1bInfo(const std::string& input_filename)
  :input_(boost::make_shared<HdfFile>(input_filename))
{
  // All work done in init    
}

std::vector<boost::shared_ptr<Level1b>> ExampleLevel1bInfo::level1b_list() {
  std::string dataset_name = "/observation_ids";
  TinyVector<int, 1> ds_shape = input_->read_shape<1>(dataset_name);
  TinyVector<int, 1> start, size;
  start = 0;
  size = ds_shape(0);
  Array<std::string, 1> obs_arr =
      input_->read_field<std::string, 1>(dataset_name, start, size);
  
  std::vector<boost::shared_ptr<Level1b>> level1b_series =
    std::vector<boost::shared_ptr<Level1b>>();
  for (auto obs_str : obs_arr) {
    boost::shared_ptr<Level1b> new_level1b =
      boost::make_shared<ExampleLevel1b>(input_, obs_str);
    level1b_series.push_back(new_level1b);
  }

  return level1b_series;
}
