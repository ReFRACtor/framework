// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%{

//--------------------------------------------------------------
/// The default conversion of a python sequence to a std::vector<T>
/// doesn't work correctly for boost::shared_ptr types. This
/// expression here does this correctly.
//--------------------------------------------------------------

// Code changed a little between SWIG 3 and SWIG 4, so have different
// versions depending on swig version.  
#if SWIG_VERSION < 0x040000  
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
	    PyObject *item;
	    typename std::vector< boost::shared_ptr<T > >::value_type temp2shared2;
	    while((item = PyIter_Next(iterator))) {
	      boost::shared_ptr<T> *itemp;
	      int newmem = 0;
	      int res = SWIG_ConvertPtrAndOwn(item, (void**) &itemp, 
			swig::type_info<boost::shared_ptr<T> >(),  0, &newmem);
	      if(!SWIG_IsOK(res)) {
		Py_DECREF(item);
		Py_DECREF(iterator);
		return SWIG_ERROR;
	      }
	      // Added mms
	      // Special handling if this is a director class. In that case, we
	      // don't own the underlying python object. See
	      // DirectorNotes.md for details.
	      Swig::Director* dp = dynamic_cast<Swig::Director*>(itemp->get());
	      if(dp) {
		temp2shared2.reset(itemp->get(), PythonRefPtrCleanup(dp->swig_get_self()));
		itemp = &temp2shared2;
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
#else
namespace swig {
  template <class SwigPySeq, class T>
  inline void
  assign(const SwigPySeq& swigpyseq, std::vector<boost::shared_ptr<T> >* seq) {
    // seq->assign(swigpyseq.begin(), swigpyseq.end()); // not used as not always implemented
    typedef typename SwigPySeq::value_type value_type;
    typename SwigPySeq::const_iterator it = swigpyseq.begin();
    for (;it != swigpyseq.end(); ++it) {
      value_type itemp = (value_type)(*it);
      // Added mms
      // Special handling if this is a director class. In that case, we
      // don't own the underlying python object. See
      // DirectorNotes.md for details.
      Swig::Director* dp = dynamic_cast<Swig::Director*>(itemp.get());
      if(dp) {
	// Diagnostic to make sure we end up here when we should
	//std::cerr << "Setting up PythonRefPtrCleanup\n";
      	itemp.reset(itemp.get(), PythonRefPtrCleanup(dp->swig_get_self()));
      }
      seq->insert(seq->end(),itemp);
    }
  }
}
#endif  
%}
