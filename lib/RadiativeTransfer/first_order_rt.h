#ifndef FO_RT_H
#define FO_RT_H

#include <boost/noncopyable.hpp>
#include <blitz/array.h>

#include "radiative_transfer_single_wn.h"
#include "rt_atmosphere.h"

/****************************************************************//**
 *******************************************************************/

namespace FullPhysics {

class FirstOrderRt : public RadiativeTransferSingleWn,
                     public Observer<RtAtmosphere>,
                     public boost::noncopyable {
public:
    FirstOrderRt(const boost::shared_ptr<RtAtmosphere>& Atm,
                 const blitz::Array<double, 1>& Sza, 
                 const blitz::Array<double, 1>& Zen, 
                 const blitz::Array<double, 1>& Azm,
                 bool do_solar = true,
                 bool do_thermal = false);
   
  //-----------------------------------------------------------------------
  /// For performance, we cache some data as we calculate it. This
  /// becomes stale when the Atmosphere is changed, so we observe atm
  /// and mark the cache when it changes. 
  //-----------------------------------------------------------------------
  void notify_update(const RtAtmosphere& atm) { /*TBD*/ }

  /// Only calculates intensity
  virtual int number_stokes() const { return 1; }

  /// I suppose first order doesn't really use any streams. Maybe this should be 1? Why is this defined as a required part of the interface anyways?
  virtual int number_stream() const { return 0; }

  virtual void print(std::ostream& Os, bool Short_form = false) const;

  virtual blitz::Array<double, 1> stokes_single_wn(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;
 
private:

};

}

#endif
