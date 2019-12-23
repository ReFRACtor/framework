#ifndef OPTICAL_PROP_LSI_H
#define OPTICAL_PROP_LSI_H

#include "optical_properties_wrt_rt.h"
#include "aerosol_optical.h"

namespace FullPhysics {

/****************************************************************//**
  Adds LSI specific funtionality for packing and unpacking optical 
  properties from an intermediate form that can be feed through the
  LSI algorithm.
  *******************************************************************/

class OpticalPropertiesLsi : public virtual OpticalPropertiesImpBase {
public:
    OpticalPropertiesLsi(const ArrayAd<double, 2>& packed_properties, double wavenumber, const boost::shared_ptr<AerosolOptical>& aerosol, int num_gas, int num_aerosol);

    /// Convert another optical properties class into a packed array of properties
    static ArrayAd<double, 2> pack(const boost::shared_ptr<OpticalPropertiesWrtRt>& source_properties);

    virtual ArrayAd<double, 2> gas_optical_depth_per_particle() const;
    virtual ArrayAd<double, 1> gas_optical_depth_per_layer() const { assert_init(); return gas_optical_depth_per_layer_; }

protected:

    void compute_aerosol_scattering(double wavenumber, const boost::shared_ptr<AerosolOptical>& aerosol);

};

}

#endif
