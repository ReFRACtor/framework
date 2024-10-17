// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "fp_logger.h"
%}

%{
#include "fp_exception.h"
#include <boost/algorithm/string.hpp>
  
namespace FullPhysics {
class PythonFpLogger : public LogImp {
public:
  PythonFpLogger(PyObject* Obj)
  {
    logger = Obj;
    Py_INCREF(logger);
  }
  virtual ~PythonFpLogger()
  {
    Py_DECREF(logger);
  }
  virtual void flush(log_level l)
  {
    std::string s = os.str();
    boost::trim(s);
    os.str("");
    switch(l) {
    case LogImp::DEBUG:
      PyObject_CallMethodObjArgs(logger, PyString_FromString("debug"),
				 PyString_FromString(s.c_str()), NULL);
      break;
    case LogImp::INFO:
      PyObject_CallMethodObjArgs(logger, PyString_FromString("info"),
				 PyString_FromString(s.c_str()), NULL);
      break;
    case LogImp::WARNING:
      PyObject_CallMethodObjArgs(logger, PyString_FromString("warning"),
				 PyString_FromString(s.c_str()), NULL);
      break;
    case LogImp::ERROR:
      PyObject_CallMethodObjArgs(logger, PyString_FromString("error"),
				 PyString_FromString(s.c_str()), NULL);
      break;
    case LogImp::FATAL:
      PyObject_CallMethodObjArgs(logger, PyString_FromString("critical"),
				 PyString_FromString(s.c_str()), NULL);
      break;
    default:
      throw Exception("Unknown log level");
    }
    if(PyErr_Occurred())
      throw std::runtime_error("Python error occurred:\n" + parse_python_exception());
  }
  virtual std::ostream* stream() {return 0;}
  static void turn_on_logger(PyObject* Obj)
  {
    FullPhysics::Logger::set_implementation(new FullPhysics::PythonFpLogger(Obj));
  }
  static void turn_off_logger()
  {
    FullPhysics::Logger::set_implementation(0);
  }
private:
  PyObject* logger;
};
}
%}
%base_import(logger)
%fp_shared_ptr(FullPhysics::FpLogger);
%fp_shared_ptr(FullPhysics::PythonFpLogger);

namespace FullPhysics {
class FpLogger : public LogImp {
public:
  FpLogger(int Verbosity_level = LogImp::INFO);
  virtual void flush(log_level l);
  %extend {
    static void turn_on_logger()
    {
      FullPhysics::Logger::set_implementation(new FullPhysics::FpLogger);
    }
    static void turn_off_logger()
    {
      FullPhysics::Logger::set_implementation(0);
    }
  }
};

// Note that there are issues with the lifetime of PythonFpLogger.
// The C++ code normally has destructors called after
// python has shut down, which means that our cleanup of the logger PyObject can
// seg fault (a race conditon, so doesn't always seg fault). We can avoid this
// by making the clean up part of the python cleanup, so it happens earlier than
// the C++ destructors.
// This can be done by manually cleaning up. This can be done
// automatically using python atexit, e.g.
// @atexit.register
// def python_fp_logger_cleanup():
//    PythonFpLogger.turn_off_logger()
  
class PythonFpLogger : public LogImp {
public:
  PythonFpLogger(PyObject* Obj);
  virtual void flush(log_level l);
  static void turn_on_logger(PyObject* Obj);
  static void turn_off_logger();
};
}

%init %{
  FullPhysics::Logger::set_implementation(new FullPhysics::FpLogger);
%}

