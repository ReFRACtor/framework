#ifndef EXAMPLE_L1B_INFO_H
#define EXAMPLE_L1B_INFO_H

#include <memory>
#include "hdf_file.h"
#include "level_1b_info.h"
#include "example_level_1b.h"
#include "example_observation_id.h"


namespace FullPhysics {

/****************************************************************//**
  This is an example L1B reader that reads a HDF formatted file
  that corresponds one-to-one with the expected interface.
*******************************************************************/

class ExampleLevel1bInfo: public Level1bInfo {
public:
    ExampleLevel1bInfo(const boost::shared_ptr<HdfFile>& input_file);
    ExampleLevel1bInfo(const std::string& input_filename);
    
    std::vector<std::shared_ptr<Level1b>> level1b_list();
    
private:
    boost::shared_ptr<HdfFile> input;

};
}
#endif
