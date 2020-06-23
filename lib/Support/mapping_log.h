#ifndef MAPPING_LOG_H
#define MAPPING_LOG_H

#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

#include "array_ad.h"
#include "mapping.h"


namespace FullPhysics {
/****************************************************************//**
  This class implements log encoding of coeffs for the retrieval view
  while using the real (linear) values for the forward model view.

  For additional information see docs for Mapping class.
*******************************************************************/
class MappingLog : public Mapping  {
public:
  //-----------------------------------------------------------------------
  /// Default Constructor.
  //-----------------------------------------------------------------------

  MappingLog() : map_name("log") {};
  virtual ~MappingLog() {}
  
  virtual ArrayAd<double, 1> fm_view(const ArrayAd<double, 1>& updated_coeff)
    const {
    blitz::Array<AutoDerivative<double>, 1> res(updated_coeff.rows());
    for(int i = 0; i < res.rows(); ++i)
      res(i) = std::exp(updated_coeff(i));
    return ArrayAd<double, 1>(res);
  }
  
  virtual ArrayAd<double, 1> retrieval_init
  (const ArrayAd<double, 1>& initial_coeff) const {
    blitz::Array<AutoDerivative<double>, 1> res(initial_coeff.rows());
    for(int i = 0; i < res.rows(); ++i)
      res(i) = std::log(initial_coeff(i));
    return ArrayAd<double, 1>(res);
  }

  virtual std::string name() const { return map_name; }

  virtual boost::shared_ptr<Mapping> clone() const
  { return boost::shared_ptr<Mapping>(new MappingLog()); }

private:
  std::string map_name;
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(MappingLog);

#endif
