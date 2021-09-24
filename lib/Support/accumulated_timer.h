#ifndef ACCUMULATED_TIMER_H
#define ACCUMULATED_TIMER_H
#include "printable.h"
#include "logger.h"
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <boost/timer/timer.hpp>

namespace FullPhysics {
class FunctionTimer;
class FunctionTimerR;

/****************************************************************//**
  This is a simple timer class that can be used to accumulate the time
  spent in multiple calls to a function. You get a function timer from
  this class, which keeps track of the time that object exists and
  adds it to the elapsed time.
*******************************************************************/
class AccumulatedTimer : public Printable<AccumulatedTimer> {
public:

//-----------------------------------------------------------------------
/// Create an AccumulatedTimer. The description is printed out when
/// you print this object.
//-----------------------------------------------------------------------

  AccumulatedTimer(const std::string& Desc) : elapsed_(0.0), desc(Desc) {}

//-----------------------------------------------------------------------
/// Total elapsed time.
//-----------------------------------------------------------------------
  
  double elapsed() const {return elapsed_; }
  
//-----------------------------------------------------------------------
/// Reset elapsed time to 0.
//-----------------------------------------------------------------------
  
  void reset_elapsed() {elapsed_ = 0.0;}

//-----------------------------------------------------------------------
/// Function timer
//-----------------------------------------------------------------------

  FunctionTimer function_timer(bool Auto_log = false) const;
  void print(std::ostream& Os) const 
  { Os << desc << " elapsed time " << elapsed(); }
private:
  mutable double elapsed_;
  std::string desc;
  friend class FunctionTimerR;
};

/****************************************************************//**
  Helper class for AccumulatedTimer
*******************************************************************/
class FunctionTimerR : boost::noncopyable {
public:
  FunctionTimerR(const AccumulatedTimer& At, bool Auto_log)  
    : at(At), auto_log(Auto_log) {}
  ~FunctionTimerR() 
  {
    // Time is in nanonseconds, but we accumulate as seconds
    at.elapsed_ += t.elapsed().user / 1e9;
    if(auto_log) {
      // There appears to be a bug, t.format() can throw a
      // std::bad_cast. This appears to be internal to boost, where
      // it casts its nanoseconds type to double using a
      // static_cast. This really should work, but it doesn't.
      // Work around this by just using elapsed user time. This
      // is pretty much what we want anyways.
      // Logger::info() << "Current: " << at.desc << " elapsed time " 
      // 		     << t.format() << "\n";
      Logger::info() << "Current: " << at.desc << " elapsed time " 
		     << t.elapsed().user / 1e9 << "\n";
      Logger::info() << "Total:   " << at << "\n";
    }
  }
private:
  const AccumulatedTimer& at;
  boost::timer::cpu_timer t;
  bool auto_log;
};

/****************************************************************//**
  Helper class for AccumulatedTimer
*******************************************************************/

class FunctionTimer {
public:
  FunctionTimer(const AccumulatedTimer& At, bool Auto_log) : 
    p(new FunctionTimerR(At, Auto_log)) {}
private:
  boost::shared_ptr<FunctionTimerR> p;
};

inline FunctionTimer AccumulatedTimer::function_timer(bool Auto_log) const 
{ return FunctionTimer(*this, Auto_log); }
}
#endif
