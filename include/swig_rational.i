//--------------------------------------------------------------
// Add support for boost rational.
//--------------------------------------------------------------
#ifndef SWIG_MODULE_ALREADY_DONE
%module(moduleimport="from ._swig_wrap import $module", directors="1", allprotected="1") foo
#endif

%{
#include <boost/rational.hpp>
%}

namespace boost {
  template<class I>  class rational {
  public:
    rational(I n);
    rational(I n, I d);
  };
}

%template(rational_int) boost::rational<int>;
