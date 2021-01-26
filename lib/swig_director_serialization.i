// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

//--------------------------------------------------------------
// Provide support for serializing python classes that implement
// director classes. This is more complicated because we have to
// handle python serialization through pickling and C++ through
// boost.
//
// Right now we don't have this in swig_rules. It isn't mature
// enough yet to worry about making this cross project. We'll have
// this in framework, work all the kinks out, and then later move this
// into swig_rules.
//--------------------------------------------------------------

%define %swig_director_serialization(NAME)
%pythoncode %{
from .serialization_support import SwigDirectorSerialization

IlsImpBase.__bases__ += (SwigDirectorSerialization, )
%}

// We may well move this into a central macro. But haven't worked out
// yet exactly what we want here, so we'll do all of this inline for now.
%{
#include "ils_imp_basePYTHON_wrap.h"
#include "fp_serialize_support.h"
  
namespace boost {
   namespace serialization {
     template<class Archive>
     void serialize(Archive & ar, SwigDirector_IlsImpBase & t, 
		    const unsigned int UNUSED(version))
     {
       ar &  boost::serialization::make_nvp(BOOST_PP_STRINGIZE(IlsImpBase), 
	    boost::serialization::base_object<FullPhysics::IlsImpBase>(t));
     }
     template<class Archive>
     void save_construct_data(Archive & ar, const SwigDirector_IlsImpBase* t, 
			      const unsigned int UNUSED(version))
     {
       PyObject* pickle_str = PyObject_CallMethodObjArgs
	 (t->swig_get_self(),
	  PyUnicode_FromString("pickle_no_cxx_save"), NULL);
       if(PyErr_Occurred()) {
	 FullPhysics::Exception e;
	 e << "Python error occurred:\n"
	   << parse_python_exception();
	 throw e;
       }
       // Need to allow null characters, so we include the size here.
       std::string python_object(PyBytes_AsString(pickle_str),
				 PyBytes_Size(pickle_str));
       ar & FP_NVP(python_object);
     }
     template<class Archive>
     void load_construct_data(Archive & ar, SwigDirector_IlsImpBase* t,
			      const unsigned int UNUSED(version))
     {
       std::string python_object;
       ar & FP_NVP(python_object);
       PyObject* mod = PyImport_ImportModule("pickle");
       // Need to allow null characters, so we include the size here.
       PyObject* lis = PyObject_CallMethodObjArgs(mod,
	 PyUnicode_FromString("loads"),
         PyBytes_FromStringAndSize(python_object.c_str(), python_object.size()),
         NULL);
       if(PyErr_Occurred()) {
	 FullPhysics::Exception e;
	 e << "Python error occurred:\n"
	   << parse_python_exception();
	 throw e;
       }
       PyObject* func = PyTuple_GetItem(lis, 0);
       if(PyErr_Occurred()) {
	 FullPhysics::Exception e;
	 e << "Python error occurred:\n"
	   << parse_python_exception();
	 throw e;
       }
       PyObject* arg = PyTuple_GetItem(lis, 1);
       if(PyErr_Occurred()) {
	 FullPhysics::Exception e;
	 e << "Python error occurred:\n"
	   << parse_python_exception();
	 throw e;
       }
       PyObject* obj = PyObject_Call(func, arg, NULL);
       if(PyErr_Occurred()) {
	 FullPhysics::Exception e;
	 e << "Python error occurred:\n"
	   << parse_python_exception();
	 throw e;
       }
       ::new(t)SwigDirector_IlsImpBase(obj);
       boost::shared_ptr<FullPhysics::IlsImpBase>* smartresult =
	 new boost::shared_ptr<FullPhysics::IlsImpBase>(t);
       SWIG_NewPointerObj(SWIG_as_voidptr(smartresult), SWIGTYPE_p_boost__shared_ptrT_FullPhysics__IlsImpBase_t, SWIG_POINTER_NEW | SWIG_POINTER_OWN);	 
     }
   }
}

BOOST_CLASS_EXPORT_KEY(SwigDirector_IlsImpBase);
BOOST_CLASS_EXPORT_IMPLEMENT(SwigDirector_IlsImpBase);

template void boost::serialization::serialize
(boost::archive::polymorphic_oarchive& ar, SwigDirector_IlsImpBase& t,
 const unsigned int version);
template void boost::serialization::serialize
 (boost::archive::polymorphic_iarchive& ar, SwigDirector_IlsImpBase& t,
  const unsigned int version);

template void boost::serialization::save_construct_data
(boost::archive::polymorphic_oarchive & ar, const SwigDirector_IlsImpBase* t,
 const unsigned int version);

template void boost::serialization::load_construct_data
(boost::archive::polymorphic_iarchive & ar, SwigDirector_IlsImpBase* t,
 const unsigned int version);
 
%}
%enddef
