#ifndef MAPPING_GAUSSIAN_H
#define MAPPING_GAUSSIAN_H

#include <blitz/array.h>

#include "array_ad.h"
#include "mapping_imp_base.h"
#include "pressure.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements Gaussian parameterization of coeffs for the
  retrieval view  while using the calculated values for the forward
  model view.

  Note that due to the scaling by total optical depth performed in
  this class, it is only suitable for Absorbers.

  For additional information see docs for MappingImpBase class.
*******************************************************************/
class MappingGaussian : public MappingImpBase  {
public:
    //-----------------------------------------------------------------------
    /// Default constructor, should call init.
    //-----------------------------------------------------------------------
    MappingGaussian() {};

    void init(bool Linear_AOD)
    {
        linear_total = Linear_AOD;
        if (linear_total) {
            map_name = "GaussianLinear";
        } else {
            map_name = "GaussianLog";
        }
    }

    //-----------------------------------------------------------------------
    /// Constructor.
    /// \param Press Pressure object used for defining pressure levels
    /// \param Linear_AOD Specifies whether the first gaussian parameter
    ///   that represents the desired total is in linear space if true,
    ///   logarithmic space if false.
    //-----------------------------------------------------------------------
    MappingGaussian(bool Linear_AOD)
        : linear_total(Linear_AOD)
    {
        if (linear_total) {
            map_name = "GaussianLinear";
        } else {
            map_name = "GaussianLog";
        }
    };

    virtual boost::shared_ptr<MappingImpBase> clone() const
    {
      return boost::shared_ptr<MappingImpBase>(new MappingGaussian(linear_total));
    }

    //-----------------------------------------------------------------------
    /// Whether this mapping uses a linear total parameter (alternative is log)
    //-----------------------------------------------------------------------
    virtual bool is_linear_total() const
    {
        return linear_total;
    }

    //-----------------------------------------------------------------------
    /// Total aerosol optical depth of the extinction values in component.
    //-----------------------------------------------------------------------
    virtual AutoDerivative<double> total_optical_depth(ArrayAd<double, 1> component,
            const boost::shared_ptr<Pressure>& updated_press) const
    {
        ArrayAd<double, 1> pressure_grid = updated_press->pressure_grid().value;
        AutoDerivative<double> tot_component = 0.0;
        for(int layer_idx = 0; layer_idx < updated_press->number_layer(); layer_idx++) {
          AutoDerivative<double> delta_press = (pressure_grid(layer_idx + 1) - pressure_grid(layer_idx)) / 2.0;
          tot_component += (delta_press * (component(layer_idx) + component(layer_idx + 1) ));
        }
        return tot_component;
    }
    //-----------------------------------------------------------------------
    /// Calculation of forward model view of coeffs with mapping applied
    //-----------------------------------------------------------------------
    virtual const ArrayAd<double, 1> fm_view(ArrayAd<double, 1> const& updated_coeff,
            const boost::shared_ptr<Pressure>& updated_press) const {
        ArrayAd<double, 1> component(updated_press->number_level(), updated_coeff.number_variable());

        AutoDerivative<double> desired_val;

        if (linear_total) {
        desired_val = updated_coeff(0);
        } else {
            desired_val = std::exp(updated_coeff(0));

        }

        // Don't let component value go lower than a minimum value. Not clear if this
        // is actually what we want to do, see ticket #2252 for this
        // suggestion. We might instead want to use a constrained solver.
        if(desired_val < min_desired) {
            desired_val = min_desired;
        }

        int ngaussians = int((updated_coeff.rows() - 1) / 2);

        ArrayAd<double, 1> pressure_grid = updated_press->pressure_grid().value;
        AutoDerivative<double> surface_press = updated_press->surface_pressure().value;

        component.resize(updated_press->number_level(), updated_coeff.number_variable());
        component.resize_number_variable(updated_coeff.number_variable());
        component.jacobian() = 0;

        for(int g_idx = 0; g_idx < ngaussians; g_idx++) {
            AutoDerivative<double> p0 = updated_coeff(g_idx * 2 + 1);
            AutoDerivative<double> sigma = updated_coeff(g_idx * 2 + 2);

            for(int lev = 0; lev < component.rows(); lev++) {
                AutoDerivative<double> p = pressure_grid(lev) / surface_press;
                AutoDerivative<double> g_eval = std::exp( -1 * (p - p0) * (p - p0) / (2 * sigma * sigma) );

                // Because its not that easy to init ArrayAd from python to 0.0
                if (g_idx == 0) {
                    component(lev) = g_eval;
                } else {
                    component(lev) = component(lev) + g_eval;
                }
            }
        }

        AutoDerivative<double> scaling_N;

        // Protect against a divide by zero
        if (total_optical_depth(component, updated_press).value() != 0.0) {
            scaling_N = desired_val / total_optical_depth(component, updated_press);
        } else {
            scaling_N = 0.0;
        }

        for(int lev = 0; lev < component.rows(); lev++) {
            component(lev) = component(lev) * scaling_N;
        }
        return component;
    };

    //-----------------------------------------------------------------------
    /// Calculation of initial retrieval view  of coeffs
    //-----------------------------------------------------------------------
    virtual const ArrayAd<double, 1> retrieval_init(ArrayAd<double, 1> const& initial_coeff) const {
        return initial_coeff;
    };

    //-----------------------------------------------------------------------
    /// Assigned mapping name
    //-----------------------------------------------------------------------
    virtual std::string name() const { return map_name; }

    virtual ~MappingGaussian() {};

private:
    std::string map_name;
    static const double min_desired;
    bool linear_total;

};
const double MappingGaussian::min_desired = 1e-9;
}

#endif
