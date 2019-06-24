#ifndef SWIG_TYPE_MAPPER_H
#define SWIG_TYPE_MAPPER_H
#include "swig_type_mapper_base.h"
#include "swig_to_python.h"

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
  virtual void* to_python(const boost::shared_ptr<GenericObject>& V) 
  {
    boost::shared_ptr<T> v2 = boost::dynamic_pointer_cast<T>(V);
    boost::shared_ptr<T> *v3 = v2 ? new boost::shared_ptr<T>(v2) : 0;
    return SWIG_NewPointerObj(SWIG_as_voidptr(v3), sinfo, 
			      SWIG_POINTER_OWN);
  }
  virtual ~SwigTypeMapper() {}
private:
  swig_type_info* sinfo;
};
}
#endif
