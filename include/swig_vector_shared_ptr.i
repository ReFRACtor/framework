// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%{

//--------------------------------------------------------------
/// The default conversion of a python sequence to a std::vector<T>
/// doesn't work correctly for boost::shared_ptr types. This
/// expression here does this correctly.
//--------------------------------------------------------------

namespace swig {
  template <class T>
  struct traits_asptr<std::vector<boost::shared_ptr<T> > >  {
    static int asptr(PyObject *obj, 
		     std::vector<boost::shared_ptr<T> > **vec) {
      // Long name, so shorten
      typedef std::vector<boost::shared_ptr<T> > vtype;
      if (obj == Py_None || SWIG_Python_GetSwigThis(obj)) {
	vtype *p;
	if (::SWIG_ConvertPtr(obj,(void**)&p,
			      swig::type_info<vtype>(),0) == SWIG_OK) {
	  if (vec) *vec = p;
	  return SWIG_OLDOBJ;
	}
      } else if (PySequence_Check(obj)) {
	try {
	  if (vec) {
	    vtype *pseq = new vtype();
	    PyObject *iterator = PyObject_GetIter(obj);
	    while(PyObject *item = PyIter_Next(iterator)) {
	      boost::shared_ptr<T> *itemp;
	      int newmem = 0;
	      int res = SWIG_ConvertPtrAndOwn(item, (void**) &itemp, 
			swig::type_info<boost::shared_ptr<T> >(),  0, &newmem);
	      if(!SWIG_IsOK(res)) {
		Py_DECREF(item);
		Py_DECREF(iterator);
		return SWIG_ERROR;
	      }
	      pseq->push_back(*itemp);
	      Py_DECREF(item);
	    }
	    Py_DECREF(iterator);
	    *vec = pseq;
	    return SWIG_NEWOBJ;
	  } else {
	    SwigPySequence_Cont<boost::shared_ptr<T> > swigpyseq(obj);
	    return swigpyseq.check() ? SWIG_OK : SWIG_ERROR;
	  }
	} catch (std::exception& e) {
	  if (vec) {
	    if (!PyErr_Occurred()) {
	      PyErr_SetString(PyExc_TypeError, e.what());
	    }
	  }
	  return SWIG_ERROR;
	}
      }
      return SWIG_ERROR;
    }
  };
}
%}
