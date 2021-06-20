#include "expandvars.h"
#include "fp_exception.h"
#include <wordexp.h>

using namespace FullPhysics;

//-----------------------------------------------------------------------
/// This is the equivalent of python "expandvars", it evaluates
/// environment variables, ~, etc. This is a thin wrapper around
/// "wordexp" which does a shell like expansion of the passed in
/// string.
///
/// This can handle $(shell command) also. We return the first word
/// found, so something like "$(ls)" will return the first value
/// listed by ls.
//-----------------------------------------------------------------------

std::string FullPhysics::expandvars(const std::string& Fin)
{
  wordexp_t r;
  int status = wordexp(Fin.c_str(), &r, 0);
  if(status != 0)
    throw Exception("expandvars failed for string '" + Fin + "'");
  std::string res(r.we_wordv[0]);
  wordfree(&r);
  return res;
}
