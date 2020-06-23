#ifndef MAPPING_OFFSET_H
#define MAPPING_OFFSET_H

#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

#include "array_ad.h"
#include "mapping.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a offset for the retrieval view while
  applying that offset to get the forward model view.

  For additional information see docs for Mapping class.
*******************************************************************/

class MappingOffset : public Mapping {
public:
  //-----------------------------------------------------------------------
  /// Default Constructor.
  //-----------------------------------------------------------------------

  MappingOffset(double Offset, blitz::Array<double, 1> Offsetee)
    : map_name("offset"), initial_offset_(Offset), offsetee_(Offsetee) {};

  virtual ~MappingOffset() {};

  double initial_offset() const { return initial_offset_; }
  const blitz::Array<double, 1>& offsetee() const { return offsetee_;}
  virtual ArrayAd<double, 1> fm_view(const ArrayAd<double, 1>& updated_coeff)
    const {
    blitz::Array<AutoDerivative<double>, 1> res(offsetee_.rows());
    for(int i = 0; i < res.rows(); ++i)
      res(i) = offsetee_(i) + updated_coeff(0);
    return ArrayAd<double,1>(res);
  }

  virtual ArrayAd<double, 1> retrieval_init
  (const ArrayAd<double, 1>& initial_coeff) const {
    blitz::Array<AutoDerivative<double>, 1> val(1);
    val(0) = initial_offset_;
    return ArrayAd<double, 1>(val);
  }

  virtual std::string name() const { return map_name; }

  virtual boost::shared_ptr<Mapping> clone() const
  {
    return boost::shared_ptr<Mapping>(new MappingOffset(initial_offset_,
							offsetee_));
  }

private:
  std::string map_name;
  double initial_offset_;
  blitz::Array<double, 1> offsetee_; // values being offset
  MappingOffset() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(MappingOffset);
#endif
