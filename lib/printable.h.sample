#ifndef GEOCAL_PRINTABLE_H
#define GEOCAL_PRINTABLE_H
#include "generic_object.h"
#include <sstream>

namespace GeoCal {

/****************************************************************//**
  This is a Mixin for classes that can be printed. 
  
  There are two different functions we want to use for printing. For
  normal C++ use, we want the usual "<<" notation, e.g., cout << foo.
  For use with languages such as Ruby that don't have a native stream
  type, we want to "print_to_string()", which returns a string that
  Ruby can directly work with.

  We implement both of these functions in terms of a third function,
  "print(ostream& Os)".

  Classes T that want to be printable in both C++ and other languages
  should derive from Printable<T>, and then define print(ostream& Os). 
*******************************************************************/

template<class T> class Printable : public virtual GenericObject {
public:
//-----------------------------------------------------------------------
/// Print to string. This is primarily useful for SWIG wrappers to
/// this class, e.g. a to_s method in ruby.
//-----------------------------------------------------------------------

  std::string print_to_string() const
  {
    // This reserve shouldn't really be necessary, but on a Mac
    // 10.4.11 using gcc 4.0.1, there is some kind of bug where we get
    // a "Double free" error when printing in Ruby. I never tracked
    // exactly where this occurred, but it was somewhere in the
    // iostream library when the buffer of os was resized. We just
    // reserve enough space up front so this isn't an issue. Since
    // this only gets called when printing, there shouldn't be much of
    // a performance hit in doing this.
    std::string buf("blah");
    buf.reserve(1000);
    std::ostringstream os(buf);
    ((T*) this)->print_wrapper(os);
    return os.str();
  }

  // Alternative to print, this works better with python. If this
  // is an empty string we fall back to print - so C++ classes work
  // without change but python can intercept this by overriding desc.
  virtual std::string desc() const { return ""; }

  void print_wrapper(std::ostream& Os) const
  {
    std::string d = desc();
    if(d != "")
      Os << d;
    else
      ((T*) this)->print(Os);
  }
  
//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

  friend std::ostream& operator<<(std::ostream& Os, const T& P)
  { P.print_wrapper(Os); return Os;}
};

}
#endif
