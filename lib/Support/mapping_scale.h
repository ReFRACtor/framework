#ifndef MAPPING_SCALE_H
#define MAPPING_SCALE_H

#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

#include "array_ad.h"
#include "mapping.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a scale factor for the retrieval view
  while applying that scale factor to get the forward model view.

  For additional information see docs for Mapping class.
*******************************************************************/
class MappingScale : public Mapping  {
public:
  //-----------------------------------------------------------------------
  /// Constructor.
  //-----------------------------------------------------------------------

  MappingScale(double Scale, blitz::Array<double, 1> Scalee)
    : map_name("scale"), initial_scale_factor_(Scale), scalee_(Scalee) {}

  virtual ~MappingScale() {}
  double initial_scale_factor() const { return initial_scale_factor_; }
  const blitz::Array<double, 1>& scalee() const { return scalee_;}

  virtual ArrayAd<double, 1> fm_view(const ArrayAd<double, 1>& updated_coeff)
    const {
    blitz::Array<AutoDerivative<double>, 1> res(scalee_.rows());
    for(int i = 0; i < res.rows(); ++i)
      res(i) = scalee_(i) * updated_coeff(0);
    return ArrayAd<double,1>(res);
  }

  virtual ArrayAd<double, 1>
  retrieval_init(const ArrayAd<double, 1>& UNUSED(updated_coeff)) const {
    blitz::Array<AutoDerivative<double>, 1> val(1);
    val(0) = initial_scale_factor_;
    return ArrayAd<double, 1>(val);
  }

  virtual std::string name() const { return map_name; }

  virtual boost::shared_ptr<Mapping> clone() const
  {
    return boost::shared_ptr<Mapping>(new MappingScale(initial_scale_factor_,
						       scalee_));
  }
private:
  std::string map_name;
  double initial_scale_factor_;
  blitz::Array<double, 1> scalee_; // values being scaled by scale
				  // factor
  MappingScale() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(MappingScale);

#endif
