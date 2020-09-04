#ifndef STATE_MAPPING_COMPOSITE_H
#define STATE_MAPPING_COMPOSITE_H

#include "state_mapping.h"
#include "array_ad.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements chaining together multiple mappings, applying
  the operations in the order they are supplied. It is up to the user
  to ensure that the chaining makes sense for the maps involved.

  For additional information see docs for StateMapping class.
*******************************************************************/
class StateMappingComposite : public StateMapping {
public:
  virtual ~StateMappingComposite() {}

  //-----------------------------------------------------------------------
  /// Default Constructor.
  //-----------------------------------------------------------------------

  StateMappingComposite(const std::vector<boost::shared_ptr<StateMapping> >& Mappings)
  : mappings(Mappings)
  {
  }

  //-----------------------------------------------------------------------
  /// Apply the mapped state operation from each mapping
  //-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> mapped_state(const ArrayAd<double, 1>& retrieval_values) const
  {
      ArrayAd<double, 1> mapped_values = retrieval_values.copy();
      BOOST_FOREACH(boost::shared_ptr<StateMapping> map, mappings) {
          mapped_values.reference( map->mapped_state(mapped_values) );
      }
      return mapped_values;
  }

  //-----------------------------------------------------------------------
  /// Apply modifications to the initial retrieval coefficients to each
  /// mapping
  //-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> retrieval_state(const ArrayAd<double, 1>& initial_values) const
  {
      ArrayAd<double, 1> init_values = initial_values.copy();
      BOOST_FOREACH(boost::shared_ptr<StateMapping> map, mappings) {
          init_values.reference( map->retrieval_state(init_values) );
      }
      return init_values;
  }

  //-----------------------------------------------------------------------
  /// Determing map name from the composed mappings
  //-----------------------------------------------------------------------
  
  virtual std::string name() const
  { 
      std::string map_name = "";
      for(int map_idx = 0; map_idx < mappings.size(); map_idx++) {
          map_name += mappings[map_idx]->name();
          if (map_idx != mappings.size()-1) {
              map_name += ", ";
          }
      }
      return map_name;
  }

  virtual boost::shared_ptr<StateMapping> clone() const
  {
    return boost::shared_ptr<StateMapping>(new StateMappingComposite(mappings));
  }
private:

  StateMappingComposite() { ; }

  std::vector<boost::shared_ptr<StateMapping> > mappings;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(StateMappingComposite);
#endif
