#ifndef STATE_MAPPING_INTERPOLATE_H
#define STATE_MAPPING_INTERPOLATE_H

#include "state_mapping.h"

#include "auto_derivative.h"
#include "array_ad.h"

#include "linear_interpolate.h"
#include "log_interpolate.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements mapping of a retrieval representation
  to a forward model representation through a interpolation from
  one pressure grid to another. The class is templated to accept the
  different types of interpolation operators.

  For additional information see docs for StateMapping class.
*******************************************************************/

template <template<class TX, class TY> class Interp>
class StateMappingInterpolate : public StateMapping {
public:
    virtual ~StateMappingInterpolate() {}

    //-----------------------------------------------------------------------
    /// Default Constructor.
    //-----------------------------------------------------------------------

    StateMappingInterpolate(const boost::shared_ptr<Pressure>& PressTo,
                            const boost::shared_ptr<Pressure>& PressFrom)
        : press_to(PressTo), press_from(PressFrom)
    {
    }

    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual ArrayAd<double, 1> fm_view(const ArrayAd<double, 1>& updated_coeff) const
    {
        blitz::Array<AutoDerivative<double>, 1> press_grid_from = press_from->pressure_grid().value.to_array();
        Unit from_unit = press_from->pressure_grid().units;

        blitz::Array<AutoDerivative<double>, 1> coeff_arr = updated_coeff.to_array();

        if (press_grid_from.rows() != coeff_arr.rows()) {
            Exception err;
            err << "Destination pressure grid size: " << press_grid_from.rows()
                << " does not match coefficient size: " << coeff_arr.rows();
            throw err;
        }

        Interp<AutoDerivative<double>, AutoDerivative<double> > interp(
            press_grid_from.begin(),
            press_grid_from.end(),
            coeff_arr.begin());

        blitz::Array<AutoDerivative<double>, 1> mapped_val(press_to->number_level());

        ArrayAd<double, 1> press_grid_to = press_to->pressure_grid().convert(from_unit).value;
        for (int lev_idx = 0; lev_idx < press_to->number_level(); lev_idx++) {
            mapped_val(lev_idx) = interp(press_grid_to(lev_idx));
        }

        return mapped_val;
    }

    //-----------------------------------------------------------------------
    /// Calculation of initial retrieval view  of coeffs with mapping applied
    //-----------------------------------------------------------------------

    virtual ArrayAd<double, 1> retrieval_init(const ArrayAd<double, 1>& initial_coeff) const
    {
        return initial_coeff;
    }

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------

    virtual std::string name() const
    {
        return Interp<AutoDerivative<double>, AutoDerivative<double> >::name() + " interpolate";
    }

    virtual boost::shared_ptr<StateMapping> clone() const
    {
        return boost::shared_ptr<StateMapping>(new StateMappingInterpolate(press_to, press_from));
    }

private:

    // For use in serialization
    StateMappingInterpolate() { ; }

    boost::shared_ptr<Pressure> press_to;
    boost::shared_ptr<Pressure> press_from;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

};

typedef StateMappingInterpolate<LinearInterpolate> StateMappingInterpolateLinearLinear;
typedef StateMappingInterpolate<LogLinearInterpolate> StateMappingInterpolateLogLinear;
typedef StateMappingInterpolate<LogLogInterpolate> StateMappingInterpolateLogLog;
typedef StateMappingInterpolate<LinearLogInterpolate> StateMappingInterpolateLinearLog;

}

FP_EXPORT_KEY(StateMappingInterpolateLinearLinear);
FP_EXPORT_KEY(StateMappingInterpolateLogLinear);
FP_EXPORT_KEY(StateMappingInterpolateLogLog);
FP_EXPORT_KEY(StateMappingInterpolateLinearLog);

#endif
