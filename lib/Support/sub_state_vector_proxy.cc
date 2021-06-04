#include "sub_state_vector_proxy.h"
#include "ostream_pad.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void SubStateVectorProxy::serialize(Archive & ar,
				       const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SubStateVectorObserver)
    & FP_NVP(proxied_observers);
}

FP_IMPLEMENT(SubStateVectorProxy);
#endif

//-----------------------------------------------------------------------
/// Registers the classes that will be proxied. Must be called from
/// subclasses constructor or else state vector will not be set up
/// correctly.
//-----------------------------------------------------------------------

void SubStateVectorProxy::initialize(const std::vector<boost::shared_ptr<SubStateVectorObserver> >& Proxied)
{
    proxied_observers = Proxied;

    int plen = 0;
    BOOST_FOREACH(boost::shared_ptr<SubStateVectorObserver> curr_obs, proxied_observers) {
        plen += curr_obs->sub_vector_size();
    }
    state_vector_observer_initialize(plen);
}

//-----------------------------------------------------------------------
/// Extracts the relevant portions of the passed arrays and passes
/// them to the proxied objects
//-----------------------------------------------------------------------

void SubStateVectorProxy::update_sub_state(const ArrayAd<double, 1>& Sv_sub, const blitz::Array<double, 2>& Cov_sub)
{
    int offset = 0;
    BOOST_FOREACH(boost::shared_ptr<SubStateVectorObserver> curr_obs, proxied_observers) {
        Range prox_range(offset, offset + curr_obs->sub_vector_size() - 1);
        curr_obs->update_sub_state(Sv_sub(prox_range), Cov_sub(prox_range, prox_range) );
        offset += curr_obs->sub_vector_size();
    }
}

//-----------------------------------------------------------------------
/// Extracts the relevant portions of the passed arrays and passes
/// them to the proxied objects
//-----------------------------------------------------------------------

void SubStateVectorProxy::state_vector_name_sub(blitz::Array<std::string, 1>& Sv_name) const
{
    int offset = 0;
    BOOST_FOREACH(boost::shared_ptr<SubStateVectorObserver> curr_obs, proxied_observers) {
        Range prox_range(offset, offset + curr_obs->sub_vector_size() - 1);
        blitz::Array<std::string, 1> prox_name(Sv_name(prox_range));
        curr_obs->state_vector_name_sub(prox_name);
        offset += curr_obs->sub_vector_size();
    }
}

//-----------------------------------------------------------------------
/// Extracts the relevant portions of the state vector values for the
/// proxied objects
//-----------------------------------------------------------------------
//
ArrayAd<double, 1> SubStateVectorProxy::sub_state_vector_values() const
{
    int values_len = 0;
    BOOST_FOREACH(boost::shared_ptr<SubStateVectorObserver> curr_obs, proxied_observers) {
        values_len += curr_obs->sub_vector_size();
    }

    ArrayAd<double, 1> combined_values(values_len, 0);

    int offset = 0;
    BOOST_FOREACH(boost::shared_ptr<SubStateVectorObserver> curr_obs, proxied_observers) {
        values_len += curr_obs->sub_vector_size();
        Range prox_range(offset, offset + curr_obs->sub_vector_size() - 1);
        
        ArrayAd<double, 1> curr_values(curr_obs->sub_state_vector_values());

        if(combined_values.is_constant() && !curr_values.is_constant()) {
            combined_values.resize_number_variable(combined_values.number_variable());
        }

        combined_values(prox_range) = curr_values;
    }

    return combined_values;

}
