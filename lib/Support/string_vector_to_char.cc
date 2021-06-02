#include "string_vector_to_char.h"
#include <sstream>
#include <iomanip>

using namespace FullPhysics;

StringVectorToChar::StringVectorToChar(const std::vector<std::string>& Components) : components(Components) {
    c_str = str_vec_to_c_str(Components);
    substrlen = max_substrlen(Components);
    num_substr = Components.size();
}
int StringVectorToChar::max_substrlen(std::vector<std::string> names) const {
    int max_sublen = 0;
    for (const auto& name : names) {
      if ((int) name.length() > max_sublen) {
	max_sublen = name.length();
      }
    }
    return max_sublen;
}

std::string StringVectorToChar::str_vec_to_c_str(std::vector<std::string> names) const {
    int c_substr_len = max_substrlen(names);
    std::stringstream c_ss;
    for (const auto& name : names) {
        c_ss << std::left << std::setw(c_substr_len) << name;
    }
    return c_ss.str();
}
