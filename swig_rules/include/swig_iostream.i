// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

// This defines some classes for making python objects look like
// ostream or istream objects.

%include "swig_iostream_inc.i"

%typemap(in) std::ostream& (boost::iostreams::filtering_ostream v) %{
    if($input != Py_None)
        v.push(python_fh_inserter($input));
    $1 = &v;
%}

%typemap(in) std::istream& (boost::iostreams::filtering_istream v) %{
    if($input != Py_None)
      // Don't hold any characters in buffer. We don't have an easy
      // way to put the buffer back into the python stream. May
      // revisit this if there is a performance reason, we are doing a
      // python call now for each read. But for now, just have a
      // simpler interface of not buffering.
      v.push(python_fh_inserter($input), 1);
    $1 = &v;
%}
