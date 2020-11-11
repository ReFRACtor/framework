// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "spurr_rt.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}
%base_import(rt_atmosphere)
%base_import(observer)
%base_import(radiative_transfer_single_wn)
%import "spectral_domain.i"
%import "optical_properties.i"
%fp_shared_ptr(FullPhysics::SpurrRt);

namespace FullPhysics {
class SpurrRt : public RadiativeTransferSingleWn,
                public Observer<RtAtmosphere> {
public:
  SpurrRt(const boost::shared_ptr<RtAtmosphere>& Atm,
          const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
          const blitz::Array<double, 1>& Sza, 
          const blitz::Array<double, 1>& Zen, 
          const blitz::Array<double, 1>& Azm,
          bool do_solar = true,
          bool do_thermal = false);

  %python_attribute(number_stokes, virtual int)
  %python_attribute(surface_type, virtual int)

   // Not sure why, but swig mangles the static variable here.
   // I think the problem is that classproperty isn't really a direct
   // thing in python. Replace with functions
   //static bool serialize_full_state;
  %extend {
    static bool serialize_full_state_get() { return FullPhysics::SpurrRt::serialize_full_state; }
    static void serialize_full_state_set(bool V) { FullPhysics::SpurrRt::serialize_full_state = V; }
  }
  
  virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain, int Spec_index) const;
  virtual ArrayAd<double, 2> stokes_and_jacobian (const SpectralDomain& Spec_domain, int Spec_index) const;

  virtual blitz::Array<double, 1> stokes_single_wn(double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn(double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  %pickle_serialization();
};
}
