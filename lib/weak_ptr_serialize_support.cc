#include "weak_ptr_serialize_support.h"
#include <set>
#include <mutex>

using namespace SWIG_MAPPER_NAMESPACE;

static std::set<const GenericObject*> ptr_referenced;
static std::recursive_mutex ptr_referenced_mutex;
static thread_local int ptr_serialize_depth = 0;

void SWIG_MAPPER_NAMESPACE::clear_ptr_serialized_reference()
{
  // Only the outermost PtrSerializeLock for this thread should reset
  // the mark state -- a nested serialize_write* call (e.g. triggered
  // while pickling a SWIG director object that itself holds other
  // framework objects) must not wipe out the still-in-progress outer
  // mark pass.
  if (ptr_serialize_depth <= 1)
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

SWIG_MAPPER_NAMESPACE::PtrSerializeLock::PtrSerializeLock()
{
  ptr_referenced_mutex.lock();
  ++ptr_serialize_depth;
}

SWIG_MAPPER_NAMESPACE::PtrSerializeLock::~PtrSerializeLock()
{
  --ptr_serialize_depth;
  ptr_referenced_mutex.unlock();
}
