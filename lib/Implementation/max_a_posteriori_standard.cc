#include <max_a_posteriori_standard.h>
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void MaxAPosterioriStandard::serialize(Archive & ar,
			const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ModelMeasureStandard)
    & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MaxAPosteriori);
}

FP_IMPLEMENT(MaxAPosterioriStandard);
#endif

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(MaxAPosterioriStandard, MaxAPosteriori)
.def(luabind::constructor< const boost::shared_ptr<ForwardModel>&,
                           const boost::shared_ptr<Observation>&, 
                           const boost::shared_ptr<StateVector>&,
                           const Array<double, 1>,
                           const Array<double, 2> >())
REGISTER_LUA_END()
#endif



//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

MaxAPosterioriStandard::MaxAPosterioriStandard
(const boost::shared_ptr<ForwardModel>& forward_model,
 const boost::shared_ptr<Observation>& observation, 
 const boost::shared_ptr<StateVector>& state_vector,
 const Array<double, 1> a_priori_params,
 const Array<double, 2> a_priori_cov)
: MaxAPosteriori(a_priori_params, a_priori_cov),
  ModelMeasureStandard(forward_model, observation, state_vector) 
{
  if(Xa.rows() != sv->observer_claimed_size())
    throw Exception("A priori state vector size and state vector size expected by the model are not equal. :( ");
}

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

MaxAPosterioriStandard::MaxAPosterioriStandard
(const std::vector<boost::shared_ptr<ForwardModel> >& forward_model,
 const std::vector<boost::shared_ptr<Observation> >& observation, 
 const boost::shared_ptr<StateVector>& state_vector,
 const Array<double, 1> a_priori_params,
 const Array<double, 2> a_priori_cov)
: MaxAPosteriori(a_priori_params, a_priori_cov),
  ModelMeasureStandard(forward_model, observation, state_vector) 
{
  if(Xa.rows() != sv->observer_claimed_size())
    throw Exception("A priori state vector size and state vector size expected by the model are not equal. :( ");
}
