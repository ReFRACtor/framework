#include "state_mapping_at_indexes.h"
#include "fp_serialize_support.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void StateMappingAtIndexes::serialize(Archive& ar, const unsigned int UNUSED(version))
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StateMapping) &
      FP_NVP(full_state) & FP_NVP_(retrieval_indexes);
}

FP_IMPLEMENT(StateMappingAtIndexes);

StateMappingAtIndexes::StateMappingAtIndexes(const blitz::Array<bool, 1>& flags)
{
    retrieval_indexes_.resize(count(flags));
    int ret_idx = 0;
    for(int flag_idx = 0; flag_idx < flags.rows(); flag_idx++) {
        if(flags(flag_idx)) {
            retrieval_indexes_(ret_idx) = flag_idx;
            ret_idx++;
        }
    }
}

ArrayAd<double, 1> StateMappingAtIndexes::mapped_state(const ArrayAd<double, 1>& updated_coeff) const
{
    if (full_state.rows() == 0) {
        Exception err;
        err << "Full state has not yet been intialized.";
        throw err;
    }

    int input_idx = 0;

    if(updated_coeff.number_variable() != full_state.number_variable()) {
        full_state.resize_number_variable(updated_coeff.number_variable());
    }

    for (int ret_idx = 0; ret_idx < retrieval_indexes_.rows(); ret_idx++) {
        full_state( retrieval_indexes_(ret_idx) ) = updated_coeff(input_idx);
        input_idx++;
    }

    return full_state;
}

ArrayAd<double, 1> StateMappingAtIndexes::retrieval_state(const ArrayAd<double, 1>& initial_values) const
{
    full_state.reference(initial_values.copy());

    ArrayAd<double, 1> retrieval_subset(retrieval_indexes_.rows(), initial_values.number_variable());
    
    int out_idx = 0;

    for (int ret_idx = 0; ret_idx < retrieval_indexes_.rows(); ret_idx++) {
        int inp_idx = retrieval_indexes_(ret_idx);

        if (inp_idx < 0 || inp_idx >= initial_values.rows()) {
            Exception err;
            err << "Input index: " << inp_idx << " from retrieval indexes list at index: " << ret_idx 
                << " out of bounds for initial values of size: " << initial_values.rows();
            throw err;
        }
        retrieval_subset(out_idx) = initial_values(inp_idx);
        out_idx++;
    }

    return retrieval_subset;
}

int StateMappingAtIndexes::initial_values_index(const int retrieval_state_index) const
{
    return retrieval_indexes_(retrieval_state_index);
}

#endif

