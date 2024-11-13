#include "rt_driver_compare.h"

#include "unit_test_support.h"

#include "spurr_brdf_types.h"
#include "planck.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

RtDriverCompare::RtDriverCompare(boost::shared_ptr<SpurrRtDriver>& baseline_driver, 
                                 boost::shared_ptr<SpurrRtDriver>& comparison_driver,
                                 bool do_debug)
: baseline_driver_(baseline_driver), comparison_driver_(comparison_driver), debug(do_debug)
{
    // Default values till changed by accessors
    refl_tol = 1e-7;

    jac_atm_rt_tol = 1e-7;
    jac_atm_fd_tol = 1e-7;

    jac_surf_rt_tol = 1e-7;
    jac_surf_fd_tol = 1e-7;
    
}

void RtDriverCompare::test_radiance(boost::shared_ptr<RtDriverScenarioGeometry>& geometry,
                                    boost::shared_ptr<RtDriverScenarioSurface>& surface,
                                    boost::shared_ptr<RtDriverScenarioAtmosphere>& atmposphere)
{

    // Convenience sizez
    int nlayer = atmposphere->number_layer();
    int nstokes = 1;

    int njac_atm = atmposphere->number_jacobian_parameter();
    int njac_surf = surface->number_jacobian_parameter();

    // Set up convenience ranges
    Range all = Range::all();
    Range rlay(0, nlayer - 1);

    Range rjac_atm(0, njac_atm - 1);
    Range rjac_surf(0, njac_surf - 1);

    // For iterating over stokes index
    firstIndex stoke_idx;

    Array<double, 1> refl_base(nstokes);
    Array<double, 1> refl_comp(nstokes);

    // Jacobian storage arrays
    blitz::Array<double, 3> jac_atm_base;
    blitz::Array<double, 2> jac_surf_param_base;
    blitz::Array<double, 1> jac_surf_temp_base;
    blitz::Array<double, 2> jac_atm_temp_base;

    blitz::Array<double, 3> jac_atm_comp;
    blitz::Array<double, 2> jac_surf_param_comp;
    blitz::Array<double, 1> jac_surf_temp_comp;
    blitz::Array<double, 2> jac_atm_temp_comp;

    // References to values in objects
    blitz::Array<double, 1> taug( atmposphere->od_gas().value() );
    blitz::Array<double, 1> taur( atmposphere->od_rayleigh().value() );
    blitz::Array<double, 1> taua( taua.value().copy() );

    blitz::Array<double, 1> od( atmposphere->total_od().value().shape() );
    blitz::Array<double, 1> ssa( atmposphere->total_od().value().shape() );
    blitz::Array<double, 3> pf( atmposphere->phase_function().value().copy() );

    // Set up atmosphere jacobians to be one for taug and the other for taur
    Array<double, 2> jac_normalization(njac_atm, nlayer);
    jac_normalization = 1;

    if (pert_atm.rows() > 0) {
        jac_normalization(0, all) = taug.value();
    }

    if (pert_atm.rows() > 1) {
        jac_normalization(1, all) = taur.value();
    }

    if (pert_atm.rows() > 2) {
        jac_normalization(2, all) = taua.value();
    }

    // Make copies of surface params as it will be modified
    // internally by the routines
    ArrayAd<double, 1> surface_params_comp;
    surface_params_comp = surface->surface_parameters();

    ArrayAd<double, 1> surface_params_base;
    surface_params_base = surface->surface_parameters();

    // Run lidort and first order to generate values for comparison
    baseline_driver->reflectance_and_jacobian_calculate(
            atmposphere->heights(),
            geometry->solar_zenith(), 
            geometry->observation_zenith(),
            geometry->relative_azimuth(),
            surface->surface_type(),
            surface->surface_parameters(),
            od, ssa, pf,
            refl_base,
            jac_atm_base, jac_surf_param_base, jac_surf_temp_base, jac_atm_temp_base,
            surface->black_body(),
            atmposphere->black_body()
    );

    comparison_driver->reflectance_and_jacobian_calculate(
            atmposphere->heights(),
            geometry->solar_zenith(), 
            geometry->observation_zenith(),
            geometry->relative_azimuth(),
            surface->surface_type(),
            surface->surface_parameters(),
            od, ssa, pf,
            refl_comp,
            jac_atm_comp, jac_surf_param_comp, jac_surf_temp_comp, jac_atm_temp_comp,
            surface->black_body(),
            atmposphere->black_body()
    );

    // Unnormalize jacobians

    jac_atm_base(rjac_atm, rlay, stoke_idx) /=  jac_normalization;
    jac_atm_comp(rjac_atm, rlay, stoke_idx) /=  jac_normalization;

    if(debug) {
        std::cerr << "refl_comp = " << refl_comp << std::endl
                  << "refl_base  = " << refl_base << std::endl;

    }

    BOOST_CHECK_MATRIX_CLOSE_TOL(refl_comp, refl_base, refl_tol);

    // Atmospheric jacobian comparison
    Array<double, 1> pert_atm( atmposphere->perturbations() );

    Array<double, 1> refl_fd(refl_base.rows());

    blitz::Array<double, 3> jac_atm_fd(njac_atm, nlayer, nstokes);
    jac_atm_fd = 0;

    for(int l_idx = 0; l_idx < nlayer; l_idx++) {
        for(int p_idx = 0; p_idx < pert_atm.rows(); p_idx++) {
            blitz::Array<double, 1> taug_pert( taug.value().copy() );
            blitz::Array<double, 1> taur_pert( taur.value().copy() );
            blitz::Array<double, 1> taua_pert( taua.value().copy() );

            blitz::Array<double, 1> od_pert( od.value().shape() );
            blitz::Array<double, 1> ssa_pert( ssa.value().shape() );
            blitz::Array<double, 3> pf_pert( pf.value().copy() );

            switch (p_idx) {
            case 0:
                taug_pert(l_idx) += pert_atm(p_idx);
                break;

            case 1:
                taur_pert(l_idx) += pert_atm(p_idx);
                break;

            case 2:
                taua_pert(l_idx) += pert_atm(p_idx);
                break;
            }

            od_pert = taur_pert + taug_pert + taua_pert;
            ssa_pert = (taur_pert + aer_prop_ssa * taua_pert) / od_pert;

            double ray_wt_pert = taur_pert(l_idx) / (taur_pert(l_idx) + aer_prop_ssa * taua_pert(l_idx));
            double aer_wt_pert = 1.0 - ray_wt_pert;

            pf_pert(all, l_idx, 0) = 0;
            pf_pert(0, l_idx, 0) = 1.0;
            pf_pert(2, l_idx, 0) = ray_wt_pert * ( (1.0 - depol) / (2.0 - depol) );

            for(int mom_idx = 1; mom_idx <= nmoms; mom_idx++) {
                pf_pert(mom_idx, l_idx, 0) = pf_pert(mom_idx, l_idx, 0) + aer_wt_pert * (2 * mom_idx + 1) * pow(aer_prop_asym, mom_idx);
            }

            refl_fd = comparison_driver->reflectance_calculate(
                    atmposphere->heights(),
                    geometry->solar_zenith(), 
                    geometry->observation_zenith(),
                    geometry->relative_azimuth(),
                    surface->surface_type(),
                    surface->surface_parameters(),
                    od_pert, ssa_pert, pf_pert,
                    surface->black_body(),
                    atmposphere->black_body()
            );

            // Driver returns array of size 1
            jac_atm_fd(p_idx, l_idx, stoke_idx) = (refl_fd(stoke_idx) - refl_base(stoke_idx)) / pert_atm(p_idx);
        }
    }

    if(debug_output) {
        for(int s_idx = 0; s_idx < num_stokes; s_idx++) {
            std::cerr << setprecision(8)
                      << "Stoke " << (s_idx+1) << ":" << std::Endl;
                      << "  jac_atm_comp = " << jac_atm_comp(rjac_atm, rlay, 0).transpose(1, 0)
                      << "  jac_atm_fd = " << jac_atm_fd(rjac_atm, rlay).transpose(1, 0)
                      << "  jac_atm_base = " << jac_atm_base(rjac_atm, rlay, 0).transpose(1, 0) << std::endl;
        }
    }

    for(int s_idx = 0; s_idx < num_stokes; s_idx++) {
        BOOST_CHECK_MATRIX_CLOSE_TOL(jac_atm_comp(rjac_atm, rlay, s_idx), jac_atm_base(rjac_atm, rlay, s_idx), jac_atm_rt_tol);
        BOOST_CHECK_MATRIX_CLOSE_TOL(jac_atm_base(rjac_atm, rlay, s_idx), jac_atm_fd(rjac_atm, rlay, s_idx), jac_atm_fd_tol);
    }

    if(blitz::any(surface_params.value() > 0.0)) {
        // Check surface jacobians against finite difference only if they are enabled
        blitz::Array<double, 1> jac_surf_param_fd( jac_surf_param_base.rows() );
        jac_surf_param_fd = 0.0;

        for(int p_idx = 0; p_idx < pert_surf.rows(); p_idx++) {
            blitz::Array<double, 1> surface_params_pert( surface->surface_parameters().rows() );
            surface_params_pert = surface->surface_parameters().value();
            surface_params_pert(p_idx) += pert_surf(p_idx);

            refl_fd = comparison_driver->reflectance_calculate(
                    atmposphere->heights(),
                    geometry->solar_zenith(), 
                    geometry->observation_zenith(),
                    geometry->relative_azimuth(),
                    surface->surface_type(),
                    surface_params_pert
                    od.value(), ssa.value(), pf.value(),
                    surface->black_body(),
                    atmposphere->black_body()
            );

            jac_surf_param_fd(p_idx) = (refl_fd(stoke_idx) - refl_base(stoke_idx)) / pert_surf(p_idx);

            // Adjust analytic jacobians to have same meaning as finite difference one
            jac_surf_param_comp(p_idx, stoke_idx) *= surface_params_comp.jacobian()(p_idx, stoke_idx);
            jac_surf_param_base(p_idx, stoke_idx) *= surface_params_base.jacobian()(p_idx, stoke_idx);
        }

        if(debug_output) {
            std::cerr << "jac_surf_param_comp = " << jac_surf_param_comp(rjac_surf, 0) << std::endl
                      << "jac_surf_param_fd = " << jac_surf_param_fd << std::endl
                      << "jac_surf_param_base = " << jac_surf_param_base(rjac_surf, 0) << std::endl;
        }

        for(int s_idx = 0; s_idx < num_stokes; s_idx++) {
            BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param_comp(rjac_surf, s_idx), jac_surf_param_base(rjac_surf, s_idx), jac_surf_rt_tol);
            BOOST_CHECK_MATRIX_CLOSE_TOL(jac_surf_param_base(rjac_surf, s_idx), jac_surf_param_fd, jac_surf_fd_tol);
        }
    }

}
