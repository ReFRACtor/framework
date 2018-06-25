#ifndef APPLY_INSTRUMENT_UNITS_H
#define APPLY_INSTRUMENT_UNITS_H

#include "instrument_correction.h"

namespace FullPhysics {

/****************************************************************//**
 An InstrumentCorrection class that applies units and scales 
 monochromatic radiances in those cases  when no other modifications 
 applies radiances.

 Useful for instance when doing thermal only calculations with no
 solar model which would normally apply units to the spectra.
*******************************************************************/

class ApplyInstrumentUnits : virtual public InstrumentCorrection,
                             public Printable<ApplyInstrumentUnits> {

public:

    ApplyInstrumentUnits(const Unit& units, const double scale_factor = 1.0) : spectral_units(units), scaling(scale_factor) {};

    virtual boost::shared_ptr<InstrumentCorrection> clone() const;
    virtual void apply_correction(const SpectralDomain& Pixel_grid, const std::vector<int>& Pixel_list, SpectralRange& Radiance) const;

    virtual void print(std::ostream& Os) const;

private:
    Unit spectral_units;
    double scaling;
};
}

#endif
