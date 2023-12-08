// Include code for mappings from std library to python
#ifndef SWIG_MODULE_ALREADY_DONE
#if SWIG_VERSION < 0x040000  
%module(directors="1", allprotected="1") foo
#else
%module(moduleimport="from ._swig_wrap import $module", directors="1", allprotected="1") foo
#endif
#endif
  
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

