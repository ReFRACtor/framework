// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%{
#include <boost/optional.hpp>
%}

// Don't think it is a problem assuming we also have blitz. If it is,
// we can separate this .i file in the future.
%import "swig_array_inc.i"

%typemap(in) boost::optional<int> (boost::optional<int> v) %{
    if($input == Py_None)
        v = boost::optional<int>();
    else
        v = boost::optional<int>(PyLong_AsLong($input));
    $1 = &v;
%}

%typemap(in) const boost::optional<int>& (boost::optional<int> v) %{
    if($input == Py_None)
        v = boost::optional<int>();
    else
        v = boost::optional<int>(PyLong_AsLong($input));
    $1 = &v;
%}
   
%typemap(out) boost::optional<int> %{
    if($1)
        $result = PyLong_FromLong(*$1);
    else
    {
        $result = Py_None;
        Py_INCREF(Py_None);
    }
%}

%typemap(out) const boost::optional<int>& %{
    if($1)
        $result = PyLong_FromLong(*$1);
    else
    {
        $result = Py_None;
        Py_INCREF(Py_None);
    }
%}

%typemap(in) boost::optional<double> (boost::optional<double> v) %{
    if($input == Py_None)
        v = boost::optional<double>();
    else
        v = boost::optional<double>(PyFloat_AsDouble($input));
    $1 = &v;
%}

%typemap(in) const boost::optional<double>& (boost::optional<double> v) %{
    if($input == Py_None)
        v = boost::optional<double>();
    else
        v = boost::optional<double>(PyFloat_AsDouble($input));
    $1 = &v;
%}
   
%typemap(out) boost::optional<double> %{
    if($1)
        $result = PyFloat_FromDouble(*$1);
    else
    {
        $result = Py_None;
        Py_INCREF(Py_None);
    }
%}

%typemap(out) const boost::optional<double>& %{
    if($1)
        $result = PyFloat_FromDouble(*$1);
    else
    {
        $result = Py_None;
        Py_INCREF(Py_None);
    }
%}

%typemap(in) const boost::optional<blitz::Range>& (boost::optional<blitz::Range> v) %{
    if($input == Py_None)
        v = boost::optional<blitz::Range>();
    else
        v = boost::optional<blitz::Range>($input);
    $1 = &v;
%}

%typemap(out) boost::optional<blitz::Range> %{
    if (&$1) {
        $result = SWIG_NewPointerObj(new blitz::Range((&$1)->get()), $descriptor(blitz::Range*), SWIG_POINTER_OWN | 0);
    } else {
        $result = Py_None;
        Py_INCREF(Py_None);
    }
%}

%typemap(out) const boost::optional<blitz::Range>& %{
    if (&$1) {
        $result = SWIG_NewPointerObj(new blitz::Range((&$1)->get()), $descriptor(blitz::Range*), SWIG_POINTER_OWN | 0);
    } else {
        $result = Py_None;
        Py_INCREF(Py_None);
    }
%}
