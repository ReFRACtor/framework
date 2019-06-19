#include "swig_type_mapper_base.h"

using namespace SWIG_MAPPER_NAMESPACE;

std::map<type_index, boost::shared_ptr<SwigTypeMapperBase> >
SwigTypeMapperBase::swig_type_map;
