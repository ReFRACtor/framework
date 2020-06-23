#include "weak_ptr_serialize_support.h"
#include <set>

using namespace SWIG_MAPPER_NAMESPACE;

static std::set<const GenericObject*> ptr_referenced;

void SWIG_MAPPER_NAMESPACE::clear_ptr_serialized_reference()
{
  ptr_referenced.clear();
}

void SWIG_MAPPER_NAMESPACE::add_ptr_serialized_reference(const GenericObject* P)
{
  ptr_referenced.insert(P);
}

bool SWIG_MAPPER_NAMESPACE::is_ptr_serialized(const GenericObject* P)
{
  return ptr_referenced.find(P) != ptr_referenced.end();
}
