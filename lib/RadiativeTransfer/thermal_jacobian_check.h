#ifndef THERMAL_JACOBIAN_CHECK_H
#define THERMAL_JACOBIAN_CHECK_H

#include "spurr_rt_driver.h"
#include "unit_test_support.h"
#include "planck.h"
#include "fe_disable_exception.h"
#include <vector>

namespace FullPhysics {

//-----------------------------------------------------------------------
/// Check the thermal Jacobians of a LIDORT/VLIDORT/FO driver against
/// finite differences for a non-isothermal nlayer atmosphere, for each
/// of the given surfaces, with an absorbing and a scattering atmosphere,
/// and (if do_downwelling) both upwelling and downwelling output.
/// Thermal emission is linear in the blackbody inputs, so the blackbody
/// finite differences are exact. Only the first Stokes element is checked.
/// pf_layer is the (moment, 1, element) phase function used for both
/// layers. surf_fd_index lists the surface parameters to finite-difference
/// (with relative step surf_eps), e.g. not the Cox-Munk shadowing flag.
//-----------------------------------------------------------------------

template<class Driver>
inline void check_thermal_jacobians_fd
(Driver& driver, const blitz::Array<double, 3>& pf_layer, int surface_type,
 const std::vector<blitz::Array<double, 1> >& surfaces,
 const std::vector<int>& surf_fd_index, bool do_downwelling = true,
 bool check_bb = true, double atm_tol = 1e-2, double surf_tol = 1e-6,
 double surf_eps = 1e-6, int nlayer = 2)
{
  using namespace blitz;
  FeDisableException disable_fp;
  Array<double, 1> heights(nlayer + 1), bb_atm(nlayer + 1);
  Array<double, 1> taug(nlayer), taur(nlayer);
  for(int i = 0; i <= nlayer; ++i) {
    heights(i) = 100.0 * (nlayer - i) / nlayer;
    bb_atm(i) = planck(568.69, 250.0 + 40.0 * i / nlayer);
  }
  for(int l = 0; l < nlayer; ++l)
    taug(l) = 0.3 * (l + 1) * 2.0 / nlayer;
  double sza = 0.001, azm = 0.0, zen = 0.0;
  double bb = planck(568.69, 290.0);

  ArrayAd<double, 3> pf(pf_layer.rows(), nlayer, pf_layer.extent(thirdDim), 1);
  for(int l = 0; l < nlayer; ++l)
    pf.value()(Range::all(), l, Range::all()) = pf_layer(Range::all(), 0, Range::all());
  pf.jacobian() = 0.0;

  for(int dn = 0; dn < (do_downwelling ? 2 : 1); ++dn)
  for(int is = 0; is < (int) surfaces.size(); ++is)
  for(double tr : {1e-8, 0.01}) {
    BOOST_TEST_CONTEXT("downwelling=" << dn << " surface=" << surfaces[is] << " taur=" << tr) {
    if(do_downwelling)
      driver.set_downwelling(dn);
    taur = tr;
    const Array<double, 1>& sp0 = surfaces[is];
    int nsurf = sp0.rows();

    // Jacobian variable 0: taug of each layer
    ArrayAd<double, 1> od(nlayer, 1), ssa(nlayer, 1);
    od.value() = taug + taur;
    ssa.value() = taur / od.value();
    od.jacobian()(Range::all(), 0) = 1.0;
    ssa.jacobian()(Range::all(), 0) = -taur / (od.value() * od.value());
    ArrayAd<double, 1> surface(nsurf, nsurf);
    surface.value() = sp0;
    surface.jacobian() = 0.0;
    for(int i = 0; i < nsurf; ++i)
      surface.jacobian()(i, i) = 1.0;

    Array<double, 1> refl, jac_surf_temp;
    Array<double, 3> jac_atm;
    Array<double, 2> jac_surf_param, jac_atm_temp;
    driver.reflectance_and_jacobian_calculate(heights, sza, azm, zen, surface_type,
        surface, od, ssa, pf, refl, jac_atm, jac_surf_param, jac_surf_temp,
        jac_atm_temp, bb, bb_atm);

    auto calc = [&](const Array<double, 1>& tg, const Array<double, 1>& sp,
                    double bbs, const Array<double, 1>& bba) {
      Array<double, 1> o(tg + taur), w(taur / (tg + taur));
      return driver.reflectance_calculate(heights, sza, azm, zen, surface_type,
                                          sp, o, w, pf.value(), bbs, bba)(0);
    };

    if(check_bb) {
      for(int lev = 0; lev < nlayer + 1; ++lev) {
        Array<double, 1> bba(bb_atm.copy());
        bba(lev) += 1e-3;
        BOOST_TEST_INFO("level " << lev);
        BOOST_CHECK_CLOSE(jac_atm_temp(lev, 0), (calc(taug, sp0, bb, bba) - refl(0)) / 1e-3, 1e-4);
      }
      BOOST_CHECK_SMALL(jac_surf_temp(0) - (calc(taug, sp0, bb + 1e-3, bb_atm) - refl(0)) / 1e-3, 1e-10);
    }
    for(int l = 0; l < nlayer; ++l) {
      Array<double, 1> tg(taug.copy());
      tg(l) += 1e-6;
      BOOST_TEST_INFO("layer " << l);
      BOOST_CHECK_CLOSE(jac_atm(0, l, 0), (calc(tg, sp0, bb, bb_atm) - refl(0)) / 1e-6, atm_tol);
    }
    // Chain the driver's surface Jacobian through the parameter transform
    // the BRDF driver applied (carried in surface.jacobian()). May be ~0
    // (e.g. downwelling without scattering), so compare absolutely.
    for(int i : surf_fd_index) {
      double an = 0;
      for(int k = 0; k < jac_surf_param.rows() && k < nsurf; ++k)
        an += jac_surf_param(k, 0) * surface.jacobian()(k, i);
      Array<double, 1> sp(sp0.copy());
      double eps = surf_eps * std::max(1.0, std::abs(sp0(i)));
      sp(i) += eps;
      BOOST_TEST_INFO("surface parameter " << i);
      BOOST_CHECK_SMALL(an - (calc(taug, sp, bb, bb_atm) - refl(0)) / eps, surf_tol);
    }
    }
  }
  if(do_downwelling)
    driver.set_downwelling(false);
}

/// Copy of pf_layer with the phase function (a1) replaced by
/// Henyey-Greenstein moments, (2l+1) g^l, for an anisotropic scatterer
inline blitz::Array<double, 3> hg_phase_function(const blitz::Array<double, 3>& pf_layer, double g = 0.7)
{
  blitz::Array<double, 3> pf(pf_layer.copy());
  for(int l = 0; l < pf.rows(); ++l)
    pf(l, 0, 0) = (2 * l + 1) * std::pow(g, l);
  return pf;
}

/// Lambertian black and reflecting surfaces
inline std::vector<blitz::Array<double, 1> > lambertian_test_surfaces()
{
  blitz::Array<double, 1> black(1), refl(1);
  black = 1e-6;
  refl = 0.3;
  return {black, refl};
}

/// Cox-Munk (non-Lambertian) surfaces: windspeed, refractive index,
/// Lambertian component, shadowing flag
inline std::vector<blitz::Array<double, 1> > coxmunk_test_surfaces()
{
  blitz::Array<double, 1> s1(4), s2(4);
  s1 = 7.0, 1.334, 0.0, 0.0;
  s2 = 3.0, 1.334, 0.05, 0.0;
  return {s1, s2};
}

}
#endif
