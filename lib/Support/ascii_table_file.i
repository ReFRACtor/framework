%include "fp_common.i"

%{
#include "ascii_table_file.h"
%}

%base_import(generic_object)

%fp_shared_ptr(FullPhysics::AsciiTableFile);

namespace FullPhysics {
class AsciiTableFile: public GenericObject {
public:
    AsciiTableFile(const std::string& filename, int skip_rows = 0);
    ~AsciiTableFile() = default;
    std::string print_to_string() const;
    %python_attribute(data, blitz::Array<double, 2>);
};
}
