// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "geocal_common.i"

%{
#include "generic_object_map.h"
%}

%base_import(generic_object)
%geocal_shared_ptr(GeoCal::GenericObjectMap);

namespace GeoCal {
class GenericObjectMap : public GenericObject {
public:
  GenericObjectMap() {}
  std::vector<std::string> keys();
  void set_generic(const std::string& Key,
		   const boost::shared_ptr<GenericObject>& V);
  boost::shared_ptr<GenericObject> get_generic(const std::string& Key) const;
  std::string print_to_string() const;
  %pickle_serialization();
  %pythoncode {
def __contains__(self, key):
  return key in self.keys()
  
def __getitem__(self, key):
  return self.get_generic(key)
      
def __setitem__(self, key, v):
  return self.set_generic(key, v)
      
def __getattr__(self, key):
  v = self.get_generic(key)
  if(v is None):
    raise AttributeError()
  return v
      }
};
}

// List of things "import *" will include
%python_export("GenericObjectMap")

