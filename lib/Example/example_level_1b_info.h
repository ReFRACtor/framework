#ifndef EXAMPLE_L1B_INFO_H
#define EXAMPLE_L1B_INFO_H

#include <memory>
#include "hdf_file.h"
#include "level_1b_info.h"
#include "example_level_1b.h"

namespace FullPhysics {

/****************************************************************//**
  This is an example L1B reader that reads a HDF formatted file
  that corresponds one-to-one with the expected interface.
*******************************************************************/

class ExampleLevel1bInfo: public Level1bInfo {
public:
  ExampleLevel1bInfo(const boost::shared_ptr<HdfFile>& input_file);
  ExampleLevel1bInfo(const std::string& input_filename);
    
  std::vector<boost::shared_ptr<Level1b>> level1b_list();
  const boost::shared_ptr<HdfFile>& input() const { return input_;}
    
private:
  boost::shared_ptr<HdfFile> input_;
  ExampleLevel1bInfo() {}
  const std::string group_name = "Level1b";
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);

};
}

FP_EXPORT_KEY(ExampleLevel1bInfo);
#endif
