#ifndef SWIG_TYPE_MAPPER_H
#define SWIG_TYPE_MAPPER_H
#include "swig_type_mapper_base.h"
#include "swig_to_python.h"

std::string parse_python_exception();

namespace SWIG_MAPPER_NAMESPACE {

/****************************************************************//**
  This is the implementation of SwigTypeMapperBase.

  Important - you can only include this header file if you have
  included the swig header stuff, e.g. normally only if you are using
  this in swig code. For geocal, this just gets included by
  geocal_shared_ptr.i. If you try to include this in normal C++ you
  will likely get an error since things like SWIG_TypeQuery aren't
  defined. 
*******************************************************************/

template<class T> class SwigTypeMapper : public SwigTypeMapperBase {
public:
  SwigTypeMapper(const char* Typename) {
    sinfo = SWIG_TypeQuery(Typename);
  }
  // We have the definition of this function in shared_ptr_type_mapper.i
  virtual bool is_python_director_check(const boost::shared_ptr<GenericObject>& V)
    const
  {
    boost::shared_ptr<Swig::Director> v2 = boost::dynamic_pointer_cast<Swig::Director>(V);
    if(v2)
      return true;
    return false;
  }

  virtual std::string get_swig_python_save_string(const boost::shared_ptr<GenericObject>& V) const
  {
    boost::shared_ptr<Swig::Director> v2 = boost::dynamic_pointer_cast<Swig::Director>(V);
    if(!v2)
      return "";
    PyObject* pobj = v2->swig_get_self();
    PyObject* this_save = PyObject_GetAttr(pobj, PyString_FromString("this"));
    PyObject_SetAttr(pobj, PyString_FromString("this"), Py_None);
    std::string python_object = cpickle_dumps(pobj);
    PyObject_SetAttr(pobj, PyString_FromString("this"), this_save);
    return python_object;
  }
  
  virtual void* to_python(const boost::shared_ptr<GenericObject>& V) const
  {
    boost::shared_ptr<T> v2 = boost::dynamic_pointer_cast<T>(V);
    boost::shared_ptr<T> *v3 = v2 ? new boost::shared_ptr<T>(v2) : 0;
    return SWIG_NewPointerObj(SWIG_as_voidptr(v3), sinfo, 
			      SWIG_POINTER_OWN);
  }
  virtual ~SwigTypeMapper() {}
private:
  swig_type_info* sinfo;
  inline PyObject* cpickle_module() const
  {
    static PyObject* mod = 0;
    if(!mod)
      mod = PyImport_ImportModule("pickle");
    return mod;
  }
  inline std::string cpickle_dumps(PyObject* obj) const
  {
    PyObject* res = PyObject_CallMethodObjArgs(cpickle_module(),
					       PyString_FromString("dumps"),
					       obj, NULL);
    if(PyErr_Occurred()) {
      throw std::runtime_error("Python error occurred:\n" + parse_python_exception());
    }
    char *buf;
    Py_ssize_t len;
    PyBytes_AsStringAndSize(res, &buf, &len);
    return std::string(buf, len);
  }
  inline PyObject* cpickle_loads(const std::string& S) const
  {
    PyObject* res = PyObject_CallMethodObjArgs(cpickle_module(),
					       PyString_FromString("loads"),
					       PyBytes_FromStringAndSize(S.c_str(), S.size()), 
					       NULL);
    if(PyErr_Occurred()) {
      throw std::runtime_error("Python error occurred:\n" + parse_python_exception());
    }
    return res;
  }
};
}
#endif
