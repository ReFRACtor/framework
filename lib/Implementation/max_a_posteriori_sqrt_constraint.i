// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "fp_common.i"

%{
#include "max_a_posteriori_sqrt_constraint.h"
%}

%base_import(max_a_posteriori)
%base_import(model_measure_standard)
%base_import(observer)

%import "forward_model.i"
%import "observation.i"
%import "state_vector.i"
%import "state_mapping.i"


%fp_shared_ptr(FullPhysics::MaxAPosterioriSqrtConstraint);
namespace FullPhysics {
  class MaxAPosterioriSqrtConstraint;
}

// Observers of MaxAPosterioriSqrtConstraint
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::MaxAPosterioriSqrtConstraint>)
%fp_shared_ptr(FullPhysics::Observer<FullPhysics::MaxAPosterioriSqrtConstraint>)
     

%template(ObservableMaxAPosterioriSqrtConstraint) FullPhysics::Observable<FullPhysics::MaxAPosterioriSqrtConstraint>;

%feature("director") FullPhysics::Observer<FullPhysics::MaxAPosterioriSqrtConstraint>;
%template(ObserverMaxAPosterioriSqrtConstraint) FullPhysics::Observer<FullPhysics::MaxAPosterioriSqrtConstraint>;

namespace FullPhysics {
  
class MaxAPosterioriSqrtConstraint : public MaxAPosteriori, public ModelMeasureStandard, public Observable<MaxAPosterioriSqrtConstraint> {
public:
  MaxAPosterioriSqrtConstraint(const boost::shared_ptr<ForwardModel>& fm,
          const boost::shared_ptr<Observation>& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> sqrt_constraint,
          const boost::shared_ptr<StateMapping>& in_map =
			       boost::make_shared<StateMappingLinear>());
  MaxAPosterioriSqrtConstraint(const std::vector<boost::shared_ptr<ForwardModel> >& fm,
 	  const std::vector<boost::shared_ptr<Observation> >& observation, 
          const boost::shared_ptr<StateVector>& state_vector,
          const blitz::Array<double, 1> a_priori_params,
          const blitz::Array<double, 2> sqrt_constraint,
          const boost::shared_ptr<StateMapping>& in_map = boost::make_shared<StateMappingLinear>());
  virtual void add_observer(Observer<MaxAPosterioriSqrtConstraint>& Obs);
  virtual void remove_observer(Observer<MaxAPosterioriSqrtConstraint>& Obs);
  void get_state(bool& OUTPUT,
		 blitz::Array<double, 1>& OUTPUT,
		 blitz::Array<double, 2>& OUTPUT,
		 blitz::Array<double, 1>& OUTPUT,
		 blitz::Array<double, 2>& OUTPUT,
		 blitz::Array<double, 2>& OUTPUT,
		 blitz::Array<double, 2>& OUTPUT) const;
  
  void set_state(const bool& msrmnt_is_const_v,
		 const blitz::Array<double, 1>& M_v,
		 const blitz::Array<double, 2>& K_v,
		 const blitz::Array<double, 1>& msrmnt_v,
		 const blitz::Array<double, 2>& msrmnt_jacobian_v,
		 const blitz::Array<double, 2>& K_x_v,
		 const blitz::Array<double, 2>& msrmnt_jacobian_x_v);
  
  %python_attribute(mapping, boost::shared_ptr<StateMapping>);
  %python_attribute_nonconst(jacobian_fm, blitz::Array<double,2>);
  %python_attribute_nonconst(measurement_jacobian_fm, blitz::Array<double,2>);
  %python_attribute_nonconst(model_measure_diff_jacobian_fm, blitz::Array<double,2>);
  %pickle_serialization();
};
}

// Extra code for handling boost serialization/python pickle of
// director classes
%fp_director_serialization(max_a_posteriori_sqrt_constraint, ObserverMaxAPosterioriSqrtConstraint)
