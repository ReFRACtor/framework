// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "fp_common.i"
%{
#include "hres_wrapper.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}

%import "optical_properties.i"

%base_import(radiative_transfer_single_wn)
%fp_shared_ptr(FullPhysics::HresWrapper);

namespace FullPhysics {

%feature("notabstract") HresWrapper;

class HresWrapper : public RadiativeTransferSingleWn {
public:
  HresWrapper(const boost::shared_ptr<RadiativeTransferSingleWn>& Rt);
  %python_attribute(number_stokes, virtual int)
  %python_attribute(number_stream, virtual int)
  virtual blitz::Array<double, 1> stokes_single_wn
  (double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn
  (double Wn, int Spec_index, const boost::shared_ptr<OpticalProperties>& Opt_prop = NULL) const;
  %python_attribute(rt, boost::shared_ptr<RadiativeTransfer>)
  %pickle_serialization();
};
}
