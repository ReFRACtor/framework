#ifndef GROUND_PIECEWISE_H
#define GROUND_PIECEWISE_H

#include "ground.h"
#include "array_with_unit.h"
#include "sub_state_vector_array.h"
#include "linear_interpolate.h"

namespace FullPhysics {

/****************************************************************//**
  This class implements a ground type implemented as a piecewise
  linear interpolation. This would be a single value that is different
  for each spectral point. Subclasses should add the state vector
  identifcation information for the type of ground information it
  represents. This class is absctract.
*******************************************************************/
class GroundPiecewise: public SubStateVectorArray<Ground> {

public:

    GroundPiecewise(const ArrayWithUnit<double, 1>& spectral_points,
                    const blitz::Array<double, 1>& point_values,
                    const blitz::Array<bool, 1>& retrieval_flag);

    virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

    virtual const AutoDerivative<double> value_at_point(const DoubleWithUnit wave_point) const;
    virtual const AutoDerivative<double> value_at_wavenumber(const double wn) const;

    virtual void update_sub_state_hook();

    virtual boost::shared_ptr<Ground> clone() const = 0;

    virtual std::string sub_state_identifier() const = 0;
    virtual std::string state_vector_name_i(int i) const = 0;

    virtual void print(std::ostream& Os) const = 0;

    virtual std::string desc() const = 0;

protected:

    // Store spectral grid as wavenumbers internally after conversion
    blitz::Array<double, 1> wavenumbers;

private:

    // Interpolation object is update for each state vector update
    boost::shared_ptr<LinearInterpolate<double, AutoDerivative<double> > > ground_interp;


};
}
#endif
