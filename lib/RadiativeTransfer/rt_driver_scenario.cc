#include "rt_driver_scenario.h"

#include "spurr_brdf_types.h"
#include "planck.h"

using namespace blitz;
using namespace FullPhysics;

/****************************************************************
 * RtDriverScenarioGeometry
 ****************************************************************/

RtDriverScenarioGeometry::RtDriverScenarioGeometry(double sza, double zen, double azm)
: sza_(sza), zen_(zen), azm_(azm)
{
}

/****************************************************************
 * RtDriverScenarioSurface
 ****************************************************************/

RtDriverScenarioSurface::RtDriverScenarioSurface(int surface_type, bool do_thermal)
{
    surface_type_ = surface_type;
    do_thermal_ = do_thermal;

    switch (surface_type) {
    case LAMBERTIAN:
		surface_params_.resize(1, 1);
		break;
	case COXMUNK:
		surface_params_.resize(4, 3);
    	break;
    case BREONVEG:
    case BREONSOIL:
		surface_params_.resize(5, 1);
        break;
    default:
        Exception e("Unhandled BRDF type index: ");
        e << surface_type;
        throw e;
    }

    // Set up jacobians
    surface_params_.jacobian() = 1.0;

    // Default perturbation is 1%
    jac_perturbation_ = Array<double, 1>(surface_params_.number_variable());
    jac_perturbation_ = 0.1;

    black_body_ = 0;
}

void RtDriverScenarioSurface::compute_simple_black_body(double wn, double temperature)
{
    black_body_ = planck(wn, temperature);
}

/****************************************************************
 * RtDriverScenarioAtmosphere
 ****************************************************************/

RtDriverScenarioAtmosphere::RtDriverScenarioAtmosphere(int nlayer, int nmoms, int nstokes, bool do_aerosol, bool do_thermal)

{
    int nparams = do_aerosol ? 3 : 2;

    heights_ = Array<double, 1>(nlayer+1);
    heights_ = 0;

    taug_ = ArrayAd<double, 1>(nlayer, nparams);
    taug_ = 0;

    taur_ = ArrayAd<double, 1>(nlayer, nparams);
    taur_ = 0;

    taua_ = ArrayAd<double, 1> (nlayer, nparams);
    taua_ = 0;

    pf_ = ArrayAd<double, 3> (nmoms+1, nlayer, nstokes, nparams);
    pf_ = 0;

    heights_ = Array<double, 1>(nlayer);
    heights_ = 0;

    // Default perturbation is 1%
    jac_perturbation_ = Array<double, 1>(nparams);
    jac_perturbation_ = 0.1;

    if (do_thermal) {
      black_body_ = Array<double, 1>(nlayer + 1);
   }

}

const RtDriverScenarioAtmosphere::ArrayAd<double, 1> total_od() const
{
    ArrayAd<double, 1> od(number_layer(), number_jacobian_parameter());

    for(int lay_idx = 0; lay_idx < number_layer(); lay_idx++) {
        od(lay_idx) = taur_(lay_idx) + taug_(lay_idx) + taua_(lay_idx);
    }

    return od;
}

const RtDriverScenarioAtmosphere::ArrayAd<double, 1> singe_scattering_albedo(double aer_prop_ssa) const
{
    ArrayAd<double, 1> ssa(number_layer(), number_jacobian_parameter());

    ArrayAd<double, 1> od(total_od());

    for(int lay_idx = 0; lay_idx < number_layer(); lay_idx++) {
        ssa(lay_idx) = (taur_(lay_idx) + aer_prop_ssa*taua_(lay_idx)) / od(lay_idx);
    }
}

void RtDriverScenarioAtmosphere::compute_simple_pf(double aer_prop_ssa, double aer_prop_asym, double depol)
{

    ArrayAd<double,1> ray_wt(number_layer(), number_jacobian_parameter());
    ArrayAd<double,1> aer_wt(number_layer(), number_jacobian_parameter())

    for(int lay_idx = 0; lay_idx < number_layer(); lay_idx++) {
        ray_wt(lay_idx) = taur_(lay_idx) / (taur_(lay_idx) + aer_prop_ssa * taua_(lay_idx));
        aer_wt(lay_idx) = 1.0 - ray_wt(lay_idx);
    }
    
    for(int lay_idx = 0; lay_idx < number_layer(); lay_idx++) {
        pf_(0, lay_idx, 0) = 1.0;
        pf_(2, lay_idx, 0) = ray_wt(lay_idx) * ( (1.0 - depol) / (2.0 - depol) );
  
        for(int mom_idx = 1; mom_idx <= nmoms; mom_idx++) {
            pf_(mom_idx, lay_idx, 0) = pf_(mom_idx, lay_idx, 0) + aer_wt(lay_idx) * (2*mom_idx+1) * pow(aer_prop_asym, mom_idx);
        }
    }
}

void RtDriverScenarioAtmosphere::compute_simple_black_body(double wn, double temperature) 
{
    back_body_ = planck(wn, temperature);
}

void RtDriverScenarioAtmosphere::compute_simple_height_grid()
{
    // Simple height grid evenly spaced
    heights(0) = 100;
    for(int hidx = 1; hidx < number_layer()+1; hidx++) {
      heights(hidx) = heights(hidx-1) - heights(0)/nlayer;
    }
}

void RtDriverScenarioAtmosphere::initialize_jacobians()
{
    Range all = Range::all();

    if (number_jacobian_parameter() >= 1) {
        taug.jacobian()(all, 0) = taug.value();
    }

    if (number_jacobian_parameter() >= 2) {
        taur.jacobian()(all, 1) = taur.value();
    }

    if (number_jacobian_parameter() >= 3) {
        taua.jacobian()(all, 2) = taua.value();
    }
}
