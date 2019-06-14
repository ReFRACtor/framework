// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%{
#include <boost/optional.hpp>
%}

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
