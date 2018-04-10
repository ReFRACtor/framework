#include <cost_func.h>

using namespace FullPhysics;


#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(CostFunc)
REGISTER_LUA_END()
#endif



const char * const CostFunc::message_str() const
{
  switch( mssg ) {
  case SOLVED: return "SOLVED";
  case ERROR: return "ERROR";
  default: return "NONE";
  }
}
