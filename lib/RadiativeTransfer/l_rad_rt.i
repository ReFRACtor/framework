// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "l_rad_rt.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}

%base_import(radiative_transfer_single_wn)
%import "spectral_bound.i"
%import "rt_atmosphere.i"
%import "array_ad.i"
%import "l_rad_driver.i"
%fp_shared_ptr(FullPhysics::LRadRt);

namespace FullPhysics {

%feature("notabstract") LRadRt;

class LRadRt : public RadiativeTransferSingleWn {
public:
    LRadRt(const boost::shared_ptr<RadiativeTransferSingleWn>& Rt,
           const SpectralBound& Spec_bound,
           const blitz::Array<double, 1>& Sza, 
           const blitz::Array<double, 1>& Zen, 
           const blitz::Array<double, 1>& Azm, 
           bool Pure_nadir,
           bool Use_first_order_scatt_calc = true,
           bool Do_second_order = false);
  
    LRadRt(const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
           const boost::shared_ptr<RtAtmosphere>& Atm,
           const SpectralBound& Spec_bound,
           const blitz::Array<double, 1>& Sza, 
           const blitz::Array<double, 1>& Zen, 
           const blitz::Array<double, 1>& Azm, 
           bool Pure_nadir,
           int Number_stokes,
           bool Do_second_order = false,
           int Number_stream = 4);

  %python_attribute(number_stokes, virtual int)
  %python_attribute(number_stream, virtual int)
  %python_attribute(surface_type, virtual int)
  virtual blitz::Array<double, 1> stokes_single_wn(double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn(double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  %python_attribute(radiative_transfer, boost::shared_ptr<RadiativeTransfer>)
  %python_attribute(l_rad_driver, boost::shared_ptr<LRadDriver>)
  ArrayAd<double, 2> interp_z_matrix(double Wn);
};
}

