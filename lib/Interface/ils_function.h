#ifndef ILS_FUNCTION_H
#define ILS_FUNCTION_H
#include "printable.h"
#include "array_ad.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This class models an Instrument Line Shape (ILS) function. This
  returns the response around a given wave number, for a given set of
  wavenumbers. 

  It is not guaranteed that the function is normalized, the calling
  class should normalize this if needed. 
*******************************************************************/

class IlsFunction : public Printable<IlsFunction> {
public:
  virtual ~IlsFunction() {}

//-----------------------------------------------------------------------
/// Return response function.
///
/// Note that is function turns out to be a bit of a bottle neck
/// because it is called so many times. Most of the time the results
/// are the same size from one call to the next, so we pass in the
/// results rather than having this be a return value like we normally
/// do. This avoids recreating the array multiple times. We resize the
/// output, so it is fine if it doesn't happen to be the final result
/// size. But much of the time we avoid and extra allocation and
/// destruction. 
///
/// \param wn_center The wave number of the center of the response
///    function
/// \param wn The wavenumbers to return response function for.
/// \param res Return the response function for each of the wn value.
//-----------------------------------------------------------------------

  virtual void ils
  (const AutoDerivative<double>& wn_center,
   const blitz::Array<double, 1>& wn, ArrayAd<double, 1>& res) const = 0;

//-----------------------------------------------------------------------
/// Descriptive name of the band.
//-----------------------------------------------------------------------

  virtual std::string band_name() const = 0;

//-----------------------------------------------------------------------
/// In general, the name used in HDF files for a particular band is
/// similar but not identical to the more human readable band_name.
/// For example, with GOSAT we use the HDF field name "weak_co2", but
/// the band name is "WC-Band". This gives the HDF name to use.
///
/// The default implementation just returns the same string as the
/// band name.
//-----------------------------------------------------------------------

  virtual std::string hdf_band_name() const { return band_name();}
  virtual void print(std::ostream& os) const { os << "IlsFunction";}
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
};
}

FP_EXPORT_KEY(IlsFunction);
#endif
