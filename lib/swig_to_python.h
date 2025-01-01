#ifndef SWIG_TO_PYTHON_H
#define SWIG_TO_PYTHON_H
#include "swig_type_mapper_base.h"
#include <vector>

//#define PTR_DEBUG 1

namespace SWIG_MAPPER_NAMESPACE {
/****************************************************************//**
  Important - you can only include this header file if you have
  included the swig header stuff, e.g. normally only if you are using
  this in swig code. For geocal, this just gets included by
  geocal_shared_ptr.i. If you try to include this in normal C++ you
  will likely get an error since things like Swig::Director aren't
  defined. 
*******************************************************************/

//-----------------------------------------------------------------------
/// Function to map from a shared point to a python object.
//-----------------------------------------------------------------------

template<typename T> inline PyObject* 
swig_to_python(const boost::shared_ptr<T>& V)
{

//-----------------------------------------------------------------------
// If pointer is Null, return None.
//-----------------------------------------------------------------------

  if(!V) {
    Py_INCREF(Py_None);
    return Py_None;
  }

//-----------------------------------------------------------------------
// If underlying object is a python object wrapped in a
// Swig::Director, return the underlying python object
//-----------------------------------------------------------------------

  Swig::Director* d = dynamic_cast<Swig::Director*>(V.get());
  if(d) {
    PyObject* p = d->swig_get_self();
    Py_INCREF(p);
    return p;
  }

//-----------------------------------------------------------------------
// See if underlying type is registered in swig_type_map. If so, return the
// underlying type
//-----------------------------------------------------------------------

  void *p = SwigTypeMapperBase::map_to_python(V);
  if(p)
    return (PyObject*) p;

//-----------------------------------------------------------------------
// Otherwise, fall back to returning the type T.
//-----------------------------------------------------------------------

  return (PyObject*) SwigTypeMapperBase::map_to_python(V, typeid(T));
}

// The lifetime of the python director object is a little different for things
// like handling serialization, because we *create* the python
// director object at the C++ level. So we don't want the extra
// Py_INCREF. I think the only place this gets used is
// serialize_function.i, but you might have other examples of this
  
  
inline PyObject* 
swig_to_python_or_none(const boost::shared_ptr<GenericObject>& V)
{
//-----------------------------------------------------------------------
// If pointer is Null, return None.
//-----------------------------------------------------------------------

  if(!V) {
    Py_INCREF(Py_None);
    return Py_None;
  }

//-----------------------------------------------------------------------
// If underlying object is a python object wrapped in a
// Swig::Director, return the underlying python object
//-----------------------------------------------------------------------

  Swig::Director* d = dynamic_cast<Swig::Director*>(V.get());
  if(d) {
    PyObject* p = d->swig_get_self();
#ifdef PTR_DEBUG    
    std::cerr << "In swig_to_python_or_none\n";
    std::cerr << "p refcnt: " << p->ob_refcnt << " address " << p->ob_type << "\n";
#endif    
    Py_INCREF(p);
    return p;
  }

//-----------------------------------------------------------------------
// See if underlying type is registered in swig_type_map. If so, return the
// underlying type
//-----------------------------------------------------------------------

  void *p = SwigTypeMapperBase::map_to_python(V);
  if(p)
    return (PyObject*) p;

//-----------------------------------------------------------------------
// Otherwise, return Py_None
//-----------------------------------------------------------------------

  Py_INCREF(Py_None);
  return Py_None;
}

template<typename T> inline PyObject* 
swig_to_python(const boost::shared_ptr<T>* V)
{ return swig_to_python(*V); }

inline PyObject* 
swig_to_python_or_none(const boost::shared_ptr<GenericObject>* V)
{ return swig_to_python_or_none(*V); }

//-----------------------------------------------------------------------
/// Function to map from a vector shared point to a python object.
//-----------------------------------------------------------------------

template<typename T> inline PyObject* 
swig_to_python(const std::vector<boost::shared_ptr<T> >& V)
{
  PyObject* res = PyList_New(V.size());
  for(int i = 0; i < V.size(); ++i)
    PyList_SetItem(res, i, swig_to_python(V[i]));
  return res;
}

template<typename T> inline PyObject* 
swig_to_python(const std::vector<std::vector<boost::shared_ptr<T> > >& V)
{
  PyObject* res = PyList_New(V.size());
  for(int i = 0; i < V.size(); ++i)
    PyList_SetItem(res, i, swig_to_python(V[i]));
  return res;
}
  
} // End namespace
#endif
