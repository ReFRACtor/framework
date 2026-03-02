// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

//--------------------------------------------------------------
// Translate exceptions into the appropriate language exception type
//--------------------------------------------------------------

%include "exception.i"

%{
  // This is defined in swig_wrap.tmpl, so it gets put into
  // swig_wrap.cc
  #include "python_exception.h"
%}

%exception {
  try {
    $action
  } catch (Swig::DirectorException &e) { 
    SWIG_fail; 
  } catch (const PythonException& e) {
    e.restore_python_exception();
    SWIG_fail; 
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%feature("director:except") {
    if ($error != NULL) {
      throw PythonException();
    }
}

