#ifndef GROUND_EMISSIVITY_PIECEWISE_H
#define GROUND_EMISSIVITY_PIECEWISE_H

#include "ground.h"
#include "array_with_unit.h"
#include "sub_state_vector_array.h"
#include "linear_interpolate.h"

namespace FullPhysics {

/****************************************************************//**
  This class implements a emissivity as a ground type as a piecewise
  linear interpolation.
*******************************************************************/
class GroundEmissivityPiecewise: public SubStateVectorArray<Ground> {

public:

    GroundEmissivityPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                              const blitz::Array<double, 1>& emissivity_values,
                              const blitz::Array<bool, 1>& retrieval_flag);

    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

    virtual const AutoDerivative<double> emissivity(const DoubleWithUnit wave_point) const;
    virtual const AutoDerivative<double> emissivity(const double wn) const;

    virtual void update_sub_state_hook();

    virtual boost::shared_ptr<Ground> clone() const;

    virtual std::string sub_state_identifier() const
    {
        return "ground/emissivity_piecewise";
    }

    virtual std::string state_vector_name_i(int i) const;

    virtual void print(std::ostream& Os) const;

    virtual std::string desc() const
    {
        return "GroundEmissivityPiecewise";
    }

private:

    // Interpolation object is update for each state vector update
    boost::shared_ptr<LinearInterpolate<double, AutoDerivative<double> > > emiss_interp;

    // Store spectral grid as wavenumbers internally after conversion
    blitz::Array<double, 1> wavenumbers;

};
}
#endif
