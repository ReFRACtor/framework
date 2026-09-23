#ifndef RT_DRIVER_COMPARE_H
#define RT_DRIVER_COMPARE_H

#include "spurr_rt_driver.h"
#include "rt_driver_scenario.h"

namespace FullPhysics {

/****************************************************************//**
 *******************************************************************/
class RtDriverCompare : public Printable<RtDriverCompare> {
public:

    RtDriverCompare(boost::shared_ptr<SpurrRtDriver>& baseline_driver, 
                    boost::shared_ptr<SpurrRtDriver>& comparison_driver,
                    bool do_debug = false);

    void test_radiance(boost::shared_ptr<RtDriverScenarioGeometry>& geometry,
                       boost::shared_ptr<RtDriverScenarioSurface>& surface,
                       boost::shared_ptr<RtDriverScenarioAtmosphere>& atmosphere);

    // Reflectance tolerance
    int reflectance_tolerance() { return refl_tol; }
    void reflectance_tolerance(double tol) { refl_tol = tol; }

    // Atmospheric jacobian tolerances
    // 1. RT to RT
    // 1. RT to finite difference
    int jac_atm_rt_tolerance() { return jac_atm_rt_tol; }
    void jac_atm_rt_tolerance(double tol) { jac_atm_rt_tol = tol; }

    int jac_atm_fd_tolerance() { return jac_atm_fd_tol; }
    void jac_atm_fd_tolerance(double tol) { jac_atm_fd_tol = tol; }

    // Surface jacobian tolerances
    // 1. RT to RT
    // 1. RT to finite difference
    int jac_surf_rt_tolerance() { return jac_surf_rt_tol; }
    void jac_surf_rt_tolerance(double tol) { jac_surf_rt_tol = tol; }

    int jac_surf_fd_tolerance() { return jac_surf_fd_tol; }
    void jac_surf_fd_tolerance(double tol) { jac_surf_fd_tol = tol; }

    virtual void print(std::ostream& Os, bool Short_form = false) const;

private:

    bool debug;

    boost::shared_ptr<SpurrRtDriver> baseline_driver_;
    boost::shared_ptr<SpurrRtDriver> comparison_driver_;

    // Tolerances
    int refl_tol;

    int jac_atm_rt_tol;
    int jac_atm_fd_tol;

    int jac_surf_rt_tol;
    int jac_surf_fd_tol;

    void check_configuration() const;
};

}

#endif
