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
  virtual void state_vector_name(const StateVector& Sv, 
			  blitz::Array<std::string, 1>& Sv_name) const;
  %python_attribute(state_vector_start_index, int);
  %python_attribute(sub_vector_size, int);
  virtual void update_sub_state(
    const ArrayAd<double, 1>& Sv_sub,
    const blitz::Array<double, 2>& Cov_sub) = 0;
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);

  %python_attribute(sub_state_identifier, std::string);
  %python_attribute(sub_state_vector_values, ArrayAd<double, 1>);

  %extend {
      // Convert sub state vector names from a blitz array to a vector
      // the same way that is done for the full state vector names in
      // the StateVector interface
      std::vector<std::string> _sub_state_vector_names() const {
          blitz::Array<std::string, 1> names_arr($self->sub_state_vector_values().rows());
          $self->state_vector_name_sub(names_arr);
          std::vector<std::string> names_vec;
          for(int i = 0; i < names_arr.extent(blitz::firstDim); ++i) {
              names_vec.push_back(names_arr(i));
          }
          return names_vec;
      }
  }

  %pythoncode {
      @property
      def sub_state_vector_names(self):
          return self._sub_state_vector_names()
  }
 
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
