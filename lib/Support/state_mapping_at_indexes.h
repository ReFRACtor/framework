#ifndef STATE_MAPPING_AT_INDEXES_H
#define STATE_MAPPING_AT_INDEXES_H

#include "state_mapping.h"
#include "array_ad.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements subsetting a full state at certain indexes.
  Updates to the full state only occur at these indexes.

  This replaces the old flag mechanisms that were built into the
  old StateVector class.

  For additional information see docs for StateMapping class.
*******************************************************************/
class StateMappingAtIndexes : public StateMapping {
public:
    virtual ~StateMappingAtIndexes() {}

    //-----------------------------------------------------------------------
    /// Create using an array of indexes into the full state of those locations
    /// that should be included in the retrieval state
    //-----------------------------------------------------------------------

    StateMappingAtIndexes(const blitz::Array<int, 1>& indexes) : retrieval_indexes(indexes) {}

    //-----------------------------------------------------------------------
    /// Create using a boolean array that indicates which indexes should be
    /// included in the retrieval state
    //-----------------------------------------------------------------------

    StateMappingAtIndexes(const blitz::Array<bool, 1>& flags);

    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual ArrayAd<double, 1> fm_view(const ArrayAd<double, 1>& updated_coeff) const;

    //-----------------------------------------------------------------------
    /// Calculation of initial retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual ArrayAd<double, 1> retrieval_init(const ArrayAd<double, 1>& initial_coeff) const;

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() const
    {
        return "linear";
    }

    virtual boost::shared_ptr<StateMapping> clone() const
    {
        return boost::shared_ptr<StateMapping>(new StateMappingAtIndexes(retrieval_indexes, full_state));
    }
private:
    // For serialization
    StateMappingAtIndexes() {}

    StateMappingAtIndexes(const blitz::Array<int, 1>& indexes, const ArrayAd<double, 1> internal_state) : retrieval_indexes(indexes), full_state(internal_state) {}

    // Original full set of cofficients before subsetting
    mutable ArrayAd<double, 1> full_state;

    // Which indexes to include in the retrieval
    blitz::Array<int, 1> retrieval_indexes;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateMappingAtIndexes);
#endif
