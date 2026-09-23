%include "fp_common.i"

%{
#include "vlidort_rt.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
#include "spurr_interface_types.h"
%}

%base_import(observer)
%base_import(rt_atmosphere)
%base_import(spurr_rt)

%import "vlidort_brdf_driver.i"
%import "vlidort_rt_driver.i"

%fp_shared_ptr(FullPhysics::VlidortRt);

namespace FullPhysics {

// Force to be not abstract
%feature("notabstract") VlidortRt;

class VlidortRt : public SpurrRt  {
public:
  VlidortRt(const boost::shared_ptr<RtAtmosphere>& Atm,
            const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
            const blitz::Array<double, 1>& Sza,
            const blitz::Array<double, 1>& Zen,
            const blitz::Array<double, 1>& Azm,
            bool Pure_nadir,
            int Number_streams,
            int Number_moments,
            bool Do_multi_scatt_only,
            int nstokes = 3,
            bool do_solar_sources = true,
            bool do_thermal_emission = false,
            bool do_thermal_scattering = true);

  %python_attribute(number_stream, virtual int)
  %python_attribute(number_moment, int)

  %python_attribute(brdf_driver, boost::shared_ptr<VLidortBrdfDriver>)
  %python_attribute(rt_driver, boost::shared_ptr<VLidortRtDriver>)

  virtual void print(std::ostream& Os, bool Short_form = false) const;
  %pickle_serialization();
};

}
