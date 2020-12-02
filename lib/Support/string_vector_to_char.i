%include "fp_common.i"

%{
#include "string_vector_to_char.h"
%}
%base_import(generic_object)

%fp_shared_ptr(FullPhysics::StringVectorToChar);

namespace FullPhysics {
/****************************************************************//**
 This class helps manage 1d C strings created from components in
 a vector<std::string>
 *******************************************************************/
class StringVectorToChar: public virtual GenericObject {
public:
    StringVectorToChar(const std::vector<std::string>& Components);
    int substrlen; ///< Length of individual strings (with padding) within larger concat'd string
    int num_substr; ///< Number of individual strings within larger concat'd string
    std::string c_str; ///< full C concat'd string with padding to max substrs fixed length
    std::vector<std::string> components; ///< Backing container of C++ strings
};
}
