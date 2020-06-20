#include "observer_serialize_support.h"
using namespace FullPhysics;

ObserverSerializedMarker* ObserverSerializedMarker::instance_ = 0;

ObserverSerializedMarker& ObserverSerializedMarker::instance()
{
  if(!instance_)
    instance_ = new ObserverSerializedMarker;
  return *instance_;
}
