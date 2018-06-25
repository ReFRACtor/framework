%include "common.i"
%{
#include "apply_instrument_units.h"
#include "sub_state_vector_array.h"
%}

%base_import(instrument_correction)
%fp_shared_ptr(FullPhysics::ApplyInstrumentUnits);

namespace FullPhysics {

%feature("notabstract") ApplyInstrumentUnits;

class ApplyInstrumentUnits : virtual public InstrumentCorrection {
public:
    ApplyInstrumentUnits(const Unit& units, const double scale_factor = 1.0) : spectral_units(units), scaling(scale_factor) {};
   
    virtual boost::shared_ptr<InstrumentCorrection> clone() const;
    virtual void apply_correction(const SpectralDomain& Pixel_grid, const std::vector<int>& Pixel_list, SpectralRange& Radiance) const;

    virtual void print(std::ostream& Os) const;
};
}
