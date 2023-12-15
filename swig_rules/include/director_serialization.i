//--------------------------------------------------------------
// Function to add the code needed to support director class
// serialization
//--------------------------------------------------------------

%define %director_serialization(NAMESPACE, INCLUDEFILE, TYPE...)
%shared_ptr(SwigDirector_ ## TYPE)
%init {
  NAMESPACE::SwigTypeMapperBase::add(typeid(SwigDirector_ ## TYPE), boost::make_shared<NAMESPACE::SwigTypeMapper< SwigDirector_ ## TYPE > > ("boost::shared_ptr< SwigDirector" "TYPE > *"));
}

%{
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
       PyObject* pobj = d->swig_get_self();
       PyObject* this_save = PyObject_GetAttr(pobj, PyString_FromString("this"));
       PyObject_SetAttr(pobj, PyString_FromString("this"), Py_None);
       std::string python_object = cpickle_dumps(pobj);
       ar & BOOST_SERIALIZATION_NVP(python_object);
       PyObject_SetAttr(pobj, PyString_FromString("this"), this_save);
       Py_DECREF(this_save);
     }
     template<class Archive>
     void load_construct_data(Archive & ar, SwigDirector_ ## TYPE* d,
 			     const unsigned int version)
     {
       std::string python_object;
       ar & BOOST_SERIALIZATION_NVP(python_object);
       PyObject* pobj = cpickle_loads(python_object);
       ::new(d)SwigDirector_ ## TYPE(pobj);
       boost::shared_ptr<NAMESPACE::TYPE> *smartresult = d ? new boost::shared_ptr<NAMESPACE::TYPE>(d) : 0;
       PyObject* thisobj = SWIG_NewPointerObj(SWIG_as_voidptr(smartresult), SWIG_TypeQuery("boost::shared_ptr<NAMESPACE::TYPE>*"), SWIG_POINTER_NEW | SWIG_POINTER_OWN);
       // Transfer ownership of SwigDirector_ ## TYPE to pobj.
       SWIG_Python_SetSwigThis(pobj, thisobj);
       Py_DECREF(thisobj);
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
