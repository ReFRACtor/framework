#include "state_mapping_at_indexes.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StateMappingAtIndexes::serialize(Archive& ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateMapping) &
      FP_NVP(full_state) & FP_NVP(retrieval_indexes);
}

FP_IMPLEMENT(StateMappingAtIndexes);

StateMappingAtIndexes::StateMappingAtIndexes(const blitz::Array<bool, 1>& flags)
{
    retrieval_indexes.resize(count(flags));
    int ret_idx = 0;
    for(int flag_idx = 0; flag_idx < flags.rows(); flag_idx++) {
        if(flags(flag_idx)) {
            retrieval_indexes(ret_idx) = flag_idx;
            ret_idx++;
        }
    }
}

ArrayAd<double, 1> StateMappingAtIndexes::fm_view(const ArrayAd<double, 1>& updated_coeff) const
{
    if (full_state.rows() == 0) {
        Exception err;
        err << "Full state has not yet been intialized.";
        throw err;
    }

    int input_idx = 0;

    for (int ret_idx = 0; ret_idx < retrieval_indexes.rows(); ret_idx++) {
        full_state( retrieval_indexes(ret_idx) ) = updated_coeff(input_idx);
        input_idx++;
    }

    return full_state;
}

//-----------------------------------------------------------------------
/// Calculation of initial retrieval view  of coeffs with mapping applied
//-----------------------------------------------------------------------

ArrayAd<double, 1> StateMappingAtIndexes::retrieval_init(const ArrayAd<double, 1>& initial_coeff) const
{
    full_state.reference(initial_coeff.copy());

    ArrayAd<double, 1> retrieval_subset(retrieval_indexes.rows(), initial_coeff.number_variable());
    
    int out_idx = 0;

    for (int ret_idx = 0; ret_idx < retrieval_indexes.rows(); ret_idx++) {
        retrieval_subset(out_idx) = initial_coeff( retrieval_indexes(ret_idx) );
        out_idx++;
    }

    return retrieval_subset;
}


#endif

