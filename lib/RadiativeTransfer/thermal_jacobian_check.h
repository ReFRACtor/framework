#ifndef THERMAL_JACOBIAN_CHECK_H
#define THERMAL_JACOBIAN_CHECK_H

#include "multiscatt_rt_driver.h"
#include "unit_test_support.h"
#include "planck.h"
#include "fe_disable_exception.h"

namespace FullPhysics {

//-----------------------------------------------------------------------
/// Check the thermal Jacobians of a LIDORT/VLIDORT driver against finite
/// differences for a non-isothermal two-layer atmosphere, over the matrix
/// {upwelling, downwelling} x {black, reflecting surface} x {absorbing,
/// scattering}. Thermal emission is linear in the blackbody inputs, so the
/// blackbody finite differences are exact. Only the first Stokes element is
/// checked. pf_layer is the (moment, 1, element) phase function used for
/// both layers.
//-----------------------------------------------------------------------

inline void check_thermal_jacobians_fd(MultiScattRtDriver& driver,
                                       const blitz::Array<double, 3>& pf_layer,
                                       int surface_type)
{
  using namespace blitz;
  FeDisableException disable_fp;
  const int nlayer = 2;
  Array<double, 1> heights(nlayer + 1);
  heights = 100, 50, 0;
  double sza = 0.001, azm = 0.0, zen = 0.0;
  double bb = planck(568.69, 290.0);
  Array<double, 1> bb_atm(nlayer + 1);
  bb_atm = planck(568.69, 250.0), planck(568.69, 270.0), planck(568.69, 290.0);
  Array<double, 1> taug(nlayer), taur(nlayer);
  taug = 0.3, 0.6;

  ArrayAd<double, 3> pf(pf_layer.rows(), nlayer, pf_layer.extent(thirdDim), 1);
  for(int l = 0; l < nlayer; ++l)
    pf.value()(Range::all(), l, Range::all()) = pf_layer(Range::all(), 0, Range::all());
  pf.jacobian() = 0.0;

  for(int dn = 0; dn < 2; ++dn)
  for(double alb : {1e-6, 0.3})
  for(double tr : {1e-8, 0.01}) {
    BOOST_TEST_CONTEXT("downwelling=" << dn << " albedo=" << alb << " taur=" << tr) {
    driver.set_downwelling(dn);
    taur = tr;

    // Jacobian variable 0: taug of each layer
    ArrayAd<double, 1> od(nlayer, 1), ssa(nlayer, 1);
    od.value() = taug + taur;
    ssa.value() = taur / od.value();
    od.jacobian()(Range::all(), 0) = 1.0;
    ssa.jacobian()(Range::all(), 0) = -taur / (od.value() * od.value());
    ArrayAd<double, 1> surface(1, 1);
    surface.value()(0) = alb;
    surface.jacobian() = 1.0;

    Array<double, 1> refl, jac_surf_temp;
    Array<double, 3> jac_atm;
    Array<double, 2> jac_surf_param, jac_atm_temp;
    driver.reflectance_and_jacobian_calculate(heights, sza, azm, zen, surface_type,
        surface, od, ssa, pf, refl, jac_atm, jac_surf_param, jac_surf_temp,
        jac_atm_temp, bb, bb_atm);

    auto calc = [&](const Array<double, 1>& tg, double a, double bbs,
                    const Array<double, 1>& bba) {
      Array<double, 1> o(tg + taur), w(taur / (tg + taur)), sp(1);
      sp = a;
      return driver.reflectance_calculate(heights, sza, azm, zen, surface_type,
                                          sp, o, w, pf.value(), bbs, bba)(0);
    };

    for(int lev = 0; lev < nlayer + 1; ++lev) {
      Array<double, 1> bba(bb_atm.copy());
      bba(lev) += 1e-3;
      BOOST_TEST_INFO("level " << lev);
      BOOST_CHECK_CLOSE(jac_atm_temp(lev, 0), (calc(taug, alb, bb, bba) - refl(0)) / 1e-3, 1e-4);
    }
    BOOST_CHECK_SMALL(jac_surf_temp(0) - (calc(taug, alb, bb + 1e-3, bb_atm) - refl(0)) / 1e-3, 1e-10);
    for(int l = 0; l < nlayer; ++l) {
      Array<double, 1> tg(taug.copy());
      tg(l) += 1e-6;
      BOOST_TEST_INFO("layer " << l);
      BOOST_CHECK_CLOSE(jac_atm(0, l, 0), (calc(tg, alb, bb, bb_atm) - refl(0)) / 1e-6, 1e-2);
    }
    // May be ~0 (downwelling, non-scattering), so compare absolutely
    BOOST_CHECK_SMALL(jac_surf_param(0, 0) * surface.jacobian()(0, 0) -
                      (calc(taug, alb + 1e-6, bb, bb_atm) - refl(0)) / 1e-6, 1e-6);
    }
  }
  driver.set_downwelling(false);
}

}
#endif
