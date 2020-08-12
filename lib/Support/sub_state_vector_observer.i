%include "fp_common.i"

%{
#include "sub_state_vector_observer.h"
%}

%base_import(state_vector_observer)

%fp_shared_ptr(FullPhysics::SubStateVectorObserver);

namespace FullPhysics {

class SubStateVectorObserver : public StateVectorObserver {
public:
  virtual ~SubStateVectorObserver();
  virtual void notify_update(const StateVector& Sv);
  virtual void mark_used(const StateVector& Sv, 
			 blitz::Array<bool, 1>& Used) const;
  virtual void state_vector_name(const StateVector& Sv, 
			  blitz::Array<std::string, 1>& Sv_name) const;
  %python_attribute(state_vector_start_index, int);
  %python_attribute(sub_vector_size, int);
  virtual void update_sub_state(
    const ArrayAd<double, 1>& Sv_sub,
    const blitz::Array<double, 2>& Cov_sub) = 0;
  virtual void mark_used_sub(blitz::Array<bool, 1>& Used) const;
  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name)
    const;
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
  %pickle_serialization();
protected:
  SubStateVectorObserver(int Plen);
  SubStateVectorObserver();
  void state_vector_observer_initialize(int Plen);
  blitz::Array<double, 1> sv_full;
  blitz::Array<double, 2> sv_cov_full;
  blitz::Array<double, 1> sv_sub;
  blitz::Array<double, 2> sv_cov_sub;
};

}
