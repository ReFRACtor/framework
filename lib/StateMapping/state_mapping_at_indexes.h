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

    StateMappingAtIndexes(const blitz::Array<int, 1>& indexes) : retrieval_indexes_(indexes) {}

    //-----------------------------------------------------------------------
    /// Create using a boolean array that indicates which indexes should be
    /// included in the retrieval state
    //-----------------------------------------------------------------------

    StateMappingAtIndexes(const blitz::Array<bool, 1>& flags);

    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual ArrayAd<double, 1> mapped_state(const ArrayAd<double, 1>& retrieval_values) const;

    //-----------------------------------------------------------------------
    /// Calculation of initial retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual ArrayAd<double, 1> retrieval_state(const ArrayAd<double, 1>& initial_values) const;

    //-----------------------------------------------------------------------
    /// Index into initial values for each retrieval state entry
    //-----------------------------------------------------------------------
    virtual int initial_values_index(const int retrieval_state_index) const;

    //-----------------------------------------------------------------------
    /// Return the indexes used for the mapping
    //-----------------------------------------------------------------------
    
    virtual const blitz::Array<int, 1> retrieval_indexes() const { return retrieval_indexes_; }

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() const
    {
        return "at_indexes";
    }

    virtual boost::shared_ptr<StateMapping> clone() const
    {
        return boost::shared_ptr<StateMapping>(new StateMappingAtIndexes(retrieval_indexes_, full_state));
    }
private:
    // For serialization
    StateMappingAtIndexes() {}

    StateMappingAtIndexes(const blitz::Array<int, 1>& indexes, const ArrayAd<double, 1> internal_state) : retrieval_indexes_(indexes), full_state(internal_state) {}

    // Original full set of cofficients before subsetting
    mutable ArrayAd<double, 1> full_state;

    // Which indexes to include in the retrieval
    blitz::Array<int, 1> retrieval_indexes_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateMappingAtIndexes);
#endif
