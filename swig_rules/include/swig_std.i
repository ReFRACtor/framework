// Include code for mappings from std library to python
%include <std_string.i>
%include <std_vector.i>
%include <typemaps.i>
%include <std_iostream.i>

%template(vector_double) std::vector<double>;
%template(vector_string) std::vector<std::string>;
%template(vector_unsigned_char) std::vector<unsigned char>;
%template(vector_short_int) std::vector<short int>;
%template(vector_int) std::vector<int>;
%template(vector_float) std::vector<float>;

