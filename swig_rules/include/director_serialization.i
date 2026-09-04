// We use to have this in the %define below, but swig 4.2 has a bug
// where it prematurely expands the macros, even though it shouldn't.
// It is confused because we have a %{ inside of a %define.
%{
#ifndef CMAKE_KLUDGE_QUOTE  
#define CMAKE_KLUDGE_Q(x) #x
#define CMAKE_KLUDGE_QUOTE(x) CMAKE_KLUDGE_Q(x)
#endif

// CMAKE unfortunately uses a different name for the wrapper file that
// the standard SWIG convention. So we set up CMAKE ot pass in a
// CMAKE_SWIG_FILE_NAMES to get the right inclusion file.  
#ifndef CMAKE_SWIG_FILE_NAMES  
  #define CMAKE_KLUDGE_INCLUDE_HEADER(x) CMAKE_KLUDGE_QUOTE(x ## _wrap.h)
#else  
  #define CMAKE_KLUDGE_INCLUDE_HEADER(x) CMAKE_KLUDGE_QUOTE(x ## PYTHON_wrap.h)
#endif  
  
%}

//--------------------------------------------------------------
// Function to add the code needed to support director class
// serialization
//--------------------------------------------------------------

%define %director_serialization(NAMESPACE, INCLUDEFILE, BNAME, TYPE...)
%shared_ptr(SwigDirector_ ## TYPE)
%init {
  NAMESPACE::SwigTypeMapperBase::add(typeid(SwigDirector_ ## TYPE), boost::make_shared<NAMESPACE::SwigTypeMapper< SwigDirector_ ## TYPE > > ("boost::shared_ptr< SwigDirector" "TYPE > *"));
}

%{
 #include CMAKE_KLUDGE_INCLUDE_HEADER(BNAME)
  
 #include "INCLUDEFILE"
 namespace boost {
   namespace serialization {
     template<class Archive>
     void serialize(Archive& ar, SwigDirector_ ## TYPE& D, const unsigned int version) {
       ar & boost::serialization::make_nvp(BOOST_PP_STRINGIZE(TYPE),
 					  boost::serialization::base_object<NAMESPACE::TYPE>(D));
     }
     template<class Archive> 
     void save_construct_data(Archive & ar,
			      const SwigDirector_ ## TYPE* d, 
			      const unsigned int version)
     {
       static PyObject* swig_this_name = PyUnicode_FromString("this");
       PyObject* pobj = d->swig_get_self();
       PyObject* this_save = PyObject_GetAttr(pobj, swig_this_name);
       PyObject_SetAttr(pobj, swig_this_name, Py_None);
       // If cpickle_dumps or the archive write throws (e.g., the python
       // object has an unpicklable attribute), we still need to restore
       // "this" before letting the exception propagate. Otherwise pobj is
       // left permanently with this=None, breaking it for the rest of its
       // lifetime.
       try {
         std::string python_object = cpickle_dumps(pobj);
         ar & BOOST_SERIALIZATION_NVP(python_object);
       } catch (...) {
         PyObject_SetAttr(pobj, swig_this_name, this_save);
         Py_DECREF(this_save);
         throw;
       }
       PyObject_SetAttr(pobj, swig_this_name, this_save);
       Py_DECREF(this_save);
     }
     template<class Archive>
     void load_construct_data(Archive & ar, SwigDirector_ ## TYPE* d,
 			     const unsigned int version)
     {
       std::string python_object;
       ar & BOOST_SERIALIZATION_NVP(python_object);
       PyObject* pobj = cpickle_loads(python_object);
       try {
         ::new(d)SwigDirector_ ## TYPE(pobj);
         boost::shared_ptr<NAMESPACE::TYPE> *smartresult = d ? new boost::shared_ptr<NAMESPACE::TYPE>(d) : 0;
         PyObject* thisobj = SWIG_NewPointerObj(SWIG_as_voidptr(smartresult), SWIG_TypeQuery("boost::shared_ptr<NAMESPACE::TYPE>*"), SWIG_POINTER_NEW | SWIG_POINTER_OWN);
         // Transfer ownership of SwigDirector_ ## TYPE to pobj.
         SWIG_Python_SetSwigThis(pobj, thisobj);
         Py_DECREF(thisobj);
       } catch (...) {
         // Construction never completed, so pobj's new reference (from
         // cpickle_loads) was never handed off to anything -- release it
         // before propagating the exception.
         Py_DECREF(pobj);
         throw;
       }
     }
   }
}
BOOST_CLASS_EXPORT_KEY(SwigDirector_ ## TYPE);
BOOST_CLASS_EXPORT_IMPLEMENT(SwigDirector_ ## TYPE);
template void boost::serialization::serialize(boost::archive::polymorphic_oarchive& ar, SwigDirector_ ## TYPE& D, const unsigned int version);
template void boost::serialization::serialize(boost::archive::polymorphic_iarchive& ar, SwigDirector_ ## TYPE& D, const unsigned int version);
template void boost::serialization::load_construct_data
  (polymorphic_iarchive & ar, SwigDirector_ ## TYPE* d, const unsigned int version);
template void boost::serialization::save_construct_data
  (polymorphic_oarchive & ar, const SwigDirector_ ## TYPE* d, const unsigned int version);
%}
%enddef     
