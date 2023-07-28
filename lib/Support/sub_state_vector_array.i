%include "fp_common.i"

%{
#include "sub_state_vector_array.h"
%}

%base_import(sub_state_vector_observer)
%base_import(state_mapping)
%base_import(state_mapping_linear)

namespace FullPhysics {

template<class Base> class SubStateVectorArray: 
    public Base,
    public SubStateVectorObserver {
public:
  SubStateVectorArray();
  void init(const blitz::Array<double, 1>& Coeff, 
            boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
  void init(double Coeff,
            boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
  virtual ~SubStateVectorArray();
  virtual std::string state_vector_name_i(int i) const;
  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) const;
  virtual void update_sub_state(const ArrayAd<double, 1>& Sv_sub, const blitz::Array<double, 2>& Cov);
  virtual void update_sub_state_hook();
  %python_attribute(coefficient, ArrayAd<double, 1>);
  %python_attribute(mapped_state, ArrayAd<double, 1>);
  %python_attribute(state_mapping, const boost::shared_ptr<StateMapping>);
  %python_attribute(statevector_covariance, blitz::Array<double, 2>);
protected:
  ArrayAd<double, 1> coeff;
  blitz::Array<double, 2> cov;
  boost::shared_ptr<StateMapping> mapping;
};

template<class Base, class ObservableBase> class SubStateVectorArray2: 
    public Base,
    public SubStateVectorObserver {
public:
  SubStateVectorArray2();
  void init(const blitz::Array<double, 1>& Coeff, 
            boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
  void init(double Coeff,
            boost::shared_ptr<StateMapping> in_map = boost::make_shared<StateMappingLinear>());
  virtual std::string state_vector_name_i(int i) const;
  virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) const;
  virtual void update_sub_state(const ArrayAd<double, 1>& Sv_sub, const blitz::Array<double, 2>& Cov);
  virtual void update_sub_state_hook();
  %python_attribute(coefficient, ArrayAd<double, 1>);
  %python_attribute(mapped_state, ArrayAd<double, 1>);
  %python_attribute(state_mapping, const boost::shared_ptr<StateMapping>);
  %python_attribute(statevector_covariance, blitz::Array<double, 2>);
protected:
  ArrayAd<double, 1> coeff;
  blitz::Array<double, 2> cov;
  boost::shared_ptr<StateMapping> mapping;
};
}

// When we create directors, at least for SWIG 2.0.4 we need to explicitly
// list every virtual function. We have lots of classes that derive from
// SubStateVectorArray, so we have this utility macro to define all those
// virtual functions.
%define %sub_state_virtual_func(TYPE)
    void add_observer_and_keep_reference(boost::shared_ptr<Observer<TYPE> >& Obs);
    void add_cache_invalidated_observer(CacheInvalidatedObserver& Obs);
    void remove_cache_invalidated_observer(CacheInvalidatedObserver& Obs);
    virtual void add_observer(Observer<TYPE>& Obs);
    virtual void remove_observer(Observer<TYPE>& Obs);
    void clear_observers();
    virtual void update_sub_state_hook();
    virtual void print(std::ostream& Os) const;
    virtual void state_vector_name(const StateVector& Sv, blitz::Array<std::string, 1>& Sv_name) const;
    virtual void notify_update(const StateVector& Observed_object);
    virtual void notify_add(StateVector& Observed_object);
    virtual void notify_remove(StateVector& Observed_object);
    virtual void update_sub_state(const ArrayAd<double, 1>& Sv_sub, const blitz::Array<double, 2>& Cov_sub);
    %python_attribute(sub_state_identifier, std::string);
    virtual std::string state_vector_name_i(int i) const;
    virtual void state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) const;
%enddef
