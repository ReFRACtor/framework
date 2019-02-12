%include "common.i"
%{
#include "example_level_1b_info.h"
%}

%base_import(level_1b_info)
%import "hdf_file.i"
%import "level_1b.i"


%fp_shared_ptr(FullPhysics::ExampleLevel1bInfo);

namespace FullPhysics {

%feature("notabstract") ExampleLevel1bInfo;

class ExampleLevel1bInfo: public Level1bInfo {
public:
    ExampleLevel1bInfo(const boost::shared_ptr<HdfFile>& input_file);
    ExampleLevel1bInfo(const std::string& input_filename);
    
    std::vector<boost::shared_ptr<Level1b>> level1b_list();
};
}
