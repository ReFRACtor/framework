#ifndef MAPPING_LINEAR_H
#define MAPPING_LINEAR_H

#include "mapping.h"
#include "array_ad.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements linear 1-to-1 mapping and just provides the
  SubStateVectorArray's coefficients as-is.

  For additional information see docs for Mapping class.
*******************************************************************/
class MappingLinear : public Mapping {
public:
  virtual ~MappingLinear() {}

  //-----------------------------------------------------------------------
  /// Default Constructor.
  //-----------------------------------------------------------------------

  MappingLinear() : map_name("linear") {}

  //-----------------------------------------------------------------------
  /// Calculation of forward model view of coeffs with mapping applied
  //-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> fm_view
  (const ArrayAd<double, 1>& updated_coeff) const
  { return updated_coeff; }

  //-----------------------------------------------------------------------
  /// Calculation of initial retrieval view  of coeffs with mapping applied
  //-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> retrieval_init
  (const ArrayAd<double, 1>& initial_coeff) const
  { return initial_coeff; }

  //-----------------------------------------------------------------------
  /// Assigned mapping name
  //-----------------------------------------------------------------------
  
  virtual std::string name() const { return map_name; }

  virtual boost::shared_ptr<Mapping> clone() const
  {
    return boost::shared_ptr<Mapping>(new MappingLinear());
  }
private:
  std::string map_name;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(MappingLinear);
#endif
