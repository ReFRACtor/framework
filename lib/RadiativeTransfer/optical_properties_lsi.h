#ifndef OPTICAL_PROP_LSI_H
#define OPTICAL_PROP_LSI_H

#include "optical_properties_wrt_input.h"

namespace FullPhysics {

/****************************************************************//**
  Adds LSI specific funtionality for packing and unpacking optical 
  properties from an intermediate form that can be feed through the
  LSI algorithm.
  *******************************************************************/

class OpticalPropertiesLsi : public virtual OpticalPropertiesWrtInput {
public:
    OpticalPropertiesLsi(const ArrayAd<double, 2>& packed_properties, double wavenumber, const boost::shared_ptr<Aerosol>& aerosol, int num_gas, int num_aerosol);

    /// Convert another optical properties class into a packed array of properties
    static ArrayAd<double, 2> pack(const boost::shared_ptr<OpticalProperties>& source_properties);

};

}

#endif
