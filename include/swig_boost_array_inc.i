// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// Stuff to include everywhere for using boost::array. We also have 
// swig_boost_array.i which includes stuff compiled in one spot.

%include "swig_array_inc.i"

%{
#include <boost/array.hpp>
#include <boost/foreach.hpp>

template<class T, int D> inline boost::array<T, D> 
  to_boost_array(PyObject* numpy)
{
  blitz::Array<T, 1> b = to_blitz_array<T, 1>(numpy);
  if(b.rows() != D)
    throw std::runtime_error("Array not expeced size");
  boost::array<T, D> res;
  for(int i = 0; i < D; ++i)
    res[i]= b(i);
  return res;
}
%}
