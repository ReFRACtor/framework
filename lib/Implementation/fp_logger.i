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
    // Ignore errors, it isn't the end of the world if we can't log
    // to the python logger
    //if(PyErr_Occurred())
      //throw std::runtime_error("Python error occurred:\n" + parse_python_exception());
  }
  virtual std::ostream* stream() {return 0;}
  // We have a life time issue, where PyFpLogger gets destroyed after most
  // of python is done, including the logger code. We take a weakref to this,
  // the logger is allowed to be destroyed at the right time and we just have
  // a weakref that raises a ReferenceError when used. Since we also ignore python
  // errors, this works fine. We have python code below for
  // turn_on_logger that sets up the weakref.
  static void _v_turn_on_logger(PyObject* Obj)
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

class PythonFpLogger : public LogImp {
public:
  PythonFpLogger(PyObject* Obj);
  virtual void flush(log_level l);
  static void _v_turn_on_logger(PyObject* Obj);
  static void turn_off_logger();
%pythoncode {
@classmethod  
def turn_on_logger(cls, obj):
    import weakref
    return cls._v_turn_on_logger(weakref.proxy(obj))
}
};
}

%init %{
  FullPhysics::Logger::set_implementation(new FullPhysics::FpLogger);
%}

