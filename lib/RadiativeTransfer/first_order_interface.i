// This file was auto-generated

%include "fp_common.i"

%{
#include "first_order_interface.h"
%}



%fp_shared_ptr(FullPhysics::Fo_Dtwpgeometry_Master);
%fp_shared_ptr(FullPhysics::Fo_Sswpgeometry_Master);
%fp_shared_ptr(FullPhysics::Fo_Scalarss_Rtcalcs_I);
%fp_shared_ptr(FullPhysics::Fo_Scalarss_Rtcalcs_Ilps);
%fp_shared_ptr(FullPhysics::Fo_Scalarss_Spherfuncs);
%fp_shared_ptr(FullPhysics::Fo_Thermal_Rtcalcs_I);
%fp_shared_ptr(FullPhysics::Fo_Thermal_Rtcalcs_Ilps);

namespace FullPhysics {



class Fo_Dtwpgeometry_Master {

public:
  Fo_Dtwpgeometry_Master(const int& maxgeoms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& ngeoms_in, const int& nlayers_in, const int& npartials_in, const int& nfine_in);
  virtual ~Fo_Dtwpgeometry_Master();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxpartials, int&)
  %python_attribute(maxfine, int&)
  %python_attribute(dtr, double&)
  %python_attribute(eradius, double&)
  %python_attribute(do_upwelling, bool&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(do_partials, bool&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(npartials, int&)
  %python_attribute(nfine, int&)
  %python_attribute(partial_layeridx, blitz::Array<int, 1>&)
  %python_attribute(heights, blitz::Array<double, 1>&)
  %python_attribute(alpha_boa, blitz::Array<double, 1>&)
  %python_attribute(partial_heights, blitz::Array<double, 1>&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(radii, blitz::Array<double, 1>&)
  %python_attribute(losw_paths, blitz::Array<double, 2>&)
  %python_attribute(alpha, blitz::Array<double, 2>&)
  %python_attribute(sina, blitz::Array<double, 2>&)
  %python_attribute(cosa, blitz::Array<double, 2>&)
  %python_attribute(radii_p, blitz::Array<double, 1>&)
  %python_attribute(losp_paths, blitz::Array<double, 2>&)
  %python_attribute(alpha_p, blitz::Array<double, 2>&)
  %python_attribute(sina_p, blitz::Array<double, 2>&)
  %python_attribute(cosa_p, blitz::Array<double, 2>&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(radiifine, blitz::Array<double, 3>&)
  %python_attribute(alphafine, blitz::Array<double, 3>&)
  %python_attribute(sinfine, blitz::Array<double, 3>&)
  %python_attribute(cosfine, blitz::Array<double, 3>&)
  %python_attribute(nfinedivs_p, blitz::Array<int, 2>&)
  %python_attribute(xfine_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_p, blitz::Array<double, 3>&)
  %python_attribute(radiifine_p, blitz::Array<double, 3>&)
  %python_attribute(alphafine_p, blitz::Array<double, 3>&)
  %python_attribute(sinfine_p, blitz::Array<double, 3>&)
  %python_attribute(cosfine_p, blitz::Array<double, 3>&)
  %python_attribute(fail, bool&)
  std::string message() const;
  std::string trace() const;
  
  void run();
};


class Fo_Sswpgeometry_Master {

public:
  Fo_Sswpgeometry_Master(const int& maxgeoms_in, const int& maxszas_in, const int& maxvzas_in, const int& maxazms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& ngeoms_in, const int& nszas_in, const int& nvzas_in, const int& nazms_in, const int& nlayers_in, const int& nfine_in, const int& npartials_in);
  virtual ~Fo_Sswpgeometry_Master();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxszas, int&)
  %python_attribute(maxvzas, int&)
  %python_attribute(maxazms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxpartials, int&)
  %python_attribute(maxfine, int&)
  %python_attribute(do_obsgeom, bool&)
  %python_attribute(do_doublet, bool&)
  %python_attribute(do_chapman, bool&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(do_partials, bool&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nszas, int&)
  %python_attribute(nvzas, int&)
  %python_attribute(nazms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(nfine, int&)
  %python_attribute(npartials, int&)
  %python_attribute(partial_layeridx, blitz::Array<int, 1>&)
  %python_attribute(dtr, double&)
  %python_attribute(pie, double&)
  %python_attribute(vsign, double&)
  %python_attribute(eradius, double&)
  %python_attribute(nv_offset, blitz::Array<int, 1>&)
  %python_attribute(na_offset, blitz::Array<int, 2>&)
  %python_attribute(nd_offset, blitz::Array<int, 1>&)
  %python_attribute(heights, blitz::Array<double, 1>&)
  %python_attribute(partial_heights, blitz::Array<double, 1>&)
  %python_attribute(obsgeom_boa, blitz::Array<double, 2>&)
  %python_attribute(alpha_boa, blitz::Array<double, 1>&)
  %python_attribute(theta_boa, blitz::Array<double, 1>&)
  %python_attribute(phi_boa, blitz::Array<double, 1>&)
  %python_attribute(donadir, blitz::Array<bool, 1>&)
  %python_attribute(raycon, blitz::Array<double, 1>&)
  %python_attribute(mu0, blitz::Array<double, 1>&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(cosscat, blitz::Array<double, 1>&)
  %python_attribute(radii, blitz::Array<double, 1>&)
  %python_attribute(losw_paths, blitz::Array<double, 2>&)
  %python_attribute(alpha, blitz::Array<double, 2>&)
  %python_attribute(sina, blitz::Array<double, 2>&)
  %python_attribute(cosa, blitz::Array<double, 2>&)
  %python_attribute(sunpaths, blitz::Array<double, 3>&)
  %python_attribute(ntraverse, blitz::Array<int, 2>&)
  %python_attribute(chapfacs, blitz::Array<double, 3>&)
  %python_attribute(theta_all, blitz::Array<double, 2>&)
  %python_attribute(radii_p, blitz::Array<double, 1>&)
  %python_attribute(losp_paths, blitz::Array<double, 2>&)
  %python_attribute(alpha_p, blitz::Array<double, 2>&)
  %python_attribute(sina_p, blitz::Array<double, 2>&)
  %python_attribute(cosa_p, blitz::Array<double, 2>&)
  %python_attribute(sunpaths_p, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_p, blitz::Array<int, 2>&)
  %python_attribute(chapfacs_p, blitz::Array<double, 3>&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(radiifine, blitz::Array<double, 3>&)
  %python_attribute(alphafine, blitz::Array<double, 3>&)
  %python_attribute(sinfine, blitz::Array<double, 3>&)
  %python_attribute(cosfine, blitz::Array<double, 3>&)
  %python_attribute(sunpathsfine, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine, blitz::Array<int, 3>&)
  %python_attribute(nfinedivs_p, blitz::Array<int, 2>&)
  %python_attribute(xfine_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_p, blitz::Array<double, 3>&)
  %python_attribute(radiifine_p, blitz::Array<double, 3>&)
  %python_attribute(alphafine_p, blitz::Array<double, 3>&)
  %python_attribute(sinfine_p, blitz::Array<double, 3>&)
  %python_attribute(cosfine_p, blitz::Array<double, 3>&)
  %python_attribute(sunpathsfine_p, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_p, blitz::Array<int, 3>&)
  %python_attribute(fail, bool&)
  std::string message() const;
  std::string trace() const;
  
  void run();
};


class Fo_Scalarss_Rtcalcs_I {

public:
  Fo_Scalarss_Rtcalcs_I(const int& maxgeoms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& maxmoments_input_in, const int& max_user_levels_in, const int& ngeoms_in, const int& nlayers_in, const int& nmoments_input_in, const int& n_user_levels_in, const int& npartials_in);
  virtual ~Fo_Scalarss_Rtcalcs_I();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxpartials, int&)
  %python_attribute(maxfine, int&)
  %python_attribute(maxmoments_input, int&)
  %python_attribute(max_user_levels, int&)
  %python_attribute(do_deltam_scaling, bool&)
  %python_attribute(do_phasfunc, bool&)
  %python_attribute(do_partials, bool&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(flux, double&)
  %python_attribute(do_sources_dn, blitz::Array<bool, 2>&)
  %python_attribute(do_sources_dn_p, blitz::Array<bool, 2>&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(nmoments_input, int&)
  %python_attribute(n_user_levels, int&)
  %python_attribute(user_levels, blitz::Array<int, 1>&)
  %python_attribute(npartials, int&)
  %python_attribute(nfinedivs_p, blitz::Array<int, 2>&)
  %python_attribute(partial_outindex, blitz::Array<int, 1>&)
  %python_attribute(partial_outflag, blitz::Array<bool, 1>&)
  %python_attribute(partial_layeridx, blitz::Array<int, 1>&)
  %python_attribute(extinction, blitz::Array<double, 1>&)
  %python_attribute(deltaus, blitz::Array<double, 1>&)
  %python_attribute(omega, blitz::Array<double, 1>&)
  %python_attribute(truncfac, blitz::Array<double, 1>&)
  %python_attribute(phasmoms, blitz::Array<double, 2>&)
  %python_attribute(phasfunc_dn, blitz::Array<double, 2>&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(legpoly_dn, blitz::Array<double, 2>&)
  %python_attribute(losw_paths, blitz::Array<double, 2>&)
  %python_attribute(losp_paths, blitz::Array<double, 2>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(sunpaths, blitz::Array<double, 3>&)
  %python_attribute(ntraverse, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine, blitz::Array<int, 3>&)
  %python_attribute(xfine_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_p, blitz::Array<double, 3>&)
  %python_attribute(sunpaths_p, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_p, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_p, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_p, blitz::Array<int, 3>&)
  %python_attribute(intensity_dn, blitz::Array<double, 2>&)
  %python_attribute(cumsource_dn, blitz::Array<double, 2>&)
  %python_attribute(lostrans_dn, blitz::Array<double, 2>&)
  %python_attribute(do_surface_leaving, bool&)
  %python_attribute(do_water_leaving, bool&)
  %python_attribute(do_sources_up, blitz::Array<bool, 2>&)
  %python_attribute(do_sources_up_p, blitz::Array<bool, 2>&)
  %python_attribute(phasfunc_up, blitz::Array<double, 2>&)
  %python_attribute(reflec, blitz::Array<double, 1>&)
  %python_attribute(slterm, blitz::Array<double, 1>&)
  %python_attribute(mu0, blitz::Array<double, 1>&)
  %python_attribute(legpoly_up, blitz::Array<double, 2>&)
  %python_attribute(intensity_up, blitz::Array<double, 2>&)
  %python_attribute(intensity_db, blitz::Array<double, 2>&)
  %python_attribute(cumsource_up, blitz::Array<double, 2>&)
  %python_attribute(cumtrans, blitz::Array<double, 2>&)
  %python_attribute(lostrans_up, blitz::Array<double, 2>&)
  %python_attribute(do_upwelling, bool&)
  %python_attribute(do_dnwelling, bool&)
  %python_attribute(xfine_up, blitz::Array<double, 3>&)
  %python_attribute(wfine_up, blitz::Array<double, 3>&)
  %python_attribute(sunpaths_up, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_up, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_up, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_up, blitz::Array<int, 3>&)
  %python_attribute(xfine_dn, blitz::Array<double, 3>&)
  %python_attribute(wfine_dn, blitz::Array<double, 3>&)
  %python_attribute(sunpaths_dn, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_dn, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_dn, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_dn, blitz::Array<int, 3>&)
  %python_attribute(xfine_up_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_up_p, blitz::Array<double, 3>&)
  %python_attribute(sunpaths_up_p, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_up_p, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_up_p, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_up_p, blitz::Array<int, 3>&)
  %python_attribute(xfine_dn_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_dn_p, blitz::Array<double, 3>&)
  %python_attribute(sunpaths_dn_p, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_dn_p, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_dn_p, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_dn_p, blitz::Array<int, 3>&)
  
  void ss_integral_i_dn();
  void ss_integral_i_up();
  void ss_integral_i_updn();
};


class Fo_Scalarss_Rtcalcs_Ilps {

public:
  Fo_Scalarss_Rtcalcs_Ilps(const int& maxgeoms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& maxmoments_input_in, const int& max_user_levels_in, const int& max_atmoswfs_in, const int& ngeoms_in, const int& nlayers_in, const int& nmoments_input_in, const int& n_user_levels_in, const int& npartials_in, const int& max_surfacewfs_in, const int& max_sleavewfs_in, const int& n_sleavewfs_in, const int& n_surfacewfs_in);
  virtual ~Fo_Scalarss_Rtcalcs_Ilps();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxpartials, int&)
  %python_attribute(maxfine, int&)
  %python_attribute(maxmoments_input, int&)
  %python_attribute(max_user_levels, int&)
  %python_attribute(max_atmoswfs, int&)
  %python_attribute(do_deltam_scaling, bool&)
  %python_attribute(do_phasfunc, bool&)
  %python_attribute(do_partials, bool&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(flux, double&)
  %python_attribute(do_sources_dn, blitz::Array<bool, 2>&)
  %python_attribute(do_sources_dn_p, blitz::Array<bool, 2>&)
  %python_attribute(do_profilewfs, bool&)
  %python_attribute(lvaryflags, blitz::Array<bool, 1>&)
  %python_attribute(lvarynums, blitz::Array<int, 1>&)
  %python_attribute(lvarymoms, blitz::Array<bool, 2>&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(nmoments_input, int&)
  %python_attribute(n_user_levels, int&)
  %python_attribute(user_levels, blitz::Array<int, 1>&)
  %python_attribute(npartials, int&)
  %python_attribute(nfinedivs_p, blitz::Array<int, 2>&)
  %python_attribute(partial_outindex, blitz::Array<int, 1>&)
  %python_attribute(partial_outflag, blitz::Array<bool, 1>&)
  %python_attribute(partial_layeridx, blitz::Array<int, 1>&)
  %python_attribute(extinction, blitz::Array<double, 1>&)
  %python_attribute(deltaus, blitz::Array<double, 1>&)
  %python_attribute(omega, blitz::Array<double, 1>&)
  %python_attribute(truncfac, blitz::Array<double, 1>&)
  %python_attribute(phasmoms, blitz::Array<double, 2>&)
  %python_attribute(phasfunc_dn, blitz::Array<double, 2>&)
  %python_attribute(l_extinction, blitz::Array<double, 2>&)
  %python_attribute(l_deltaus, blitz::Array<double, 2>&)
  %python_attribute(l_omega, blitz::Array<double, 2>&)
  %python_attribute(l_truncfac, blitz::Array<double, 2>&)
  %python_attribute(l_phasmoms, blitz::Array<double, 3>&)
  %python_attribute(l_phasfunc_dn, blitz::Array<double, 3>&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(legpoly_dn, blitz::Array<double, 2>&)
  %python_attribute(losw_paths, blitz::Array<double, 2>&)
  %python_attribute(losp_paths, blitz::Array<double, 2>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(sunpaths, blitz::Array<double, 3>&)
  %python_attribute(ntraverse, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine, blitz::Array<int, 3>&)
  %python_attribute(xfine_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_p, blitz::Array<double, 3>&)
  %python_attribute(sunpaths_p, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_p, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_p, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_p, blitz::Array<int, 3>&)
  %python_attribute(intensity_dn, blitz::Array<double, 2>&)
  %python_attribute(lp_jacobians_dn, blitz::Array<double, 4>&)
  %python_attribute(lostrans_dn, blitz::Array<double, 2>&)
  %python_attribute(lp_lostrans_dn, blitz::Array<double, 3>&)
  %python_attribute(max_surfacewfs, int&)
  %python_attribute(max_sleavewfs, int&)
  %python_attribute(do_surface_leaving, bool&)
  %python_attribute(do_water_leaving, bool&)
  %python_attribute(do_sources_up, blitz::Array<bool, 2>&)
  %python_attribute(do_sources_up_p, blitz::Array<bool, 2>&)
  %python_attribute(do_surfacewfs, bool&)
  %python_attribute(do_sleavewfs, bool&)
  %python_attribute(n_reflecwfs, int&)
  %python_attribute(n_sleavewfs, int&)
  %python_attribute(n_surfacewfs, int&)
  %python_attribute(phasfunc_up, blitz::Array<double, 2>&)
  %python_attribute(reflec, blitz::Array<double, 1>&)
  %python_attribute(slterm, blitz::Array<double, 1>&)
  %python_attribute(l_phasfunc_up, blitz::Array<double, 3>&)
  %python_attribute(ls_reflec, blitz::Array<double, 2>&)
  %python_attribute(lssl_slterm, blitz::Array<double, 2>&)
  %python_attribute(mu0, blitz::Array<double, 1>&)
  %python_attribute(legpoly_up, blitz::Array<double, 2>&)
  %python_attribute(intensity_up, blitz::Array<double, 2>&)
  %python_attribute(intensity_db, blitz::Array<double, 2>&)
  %python_attribute(lp_jacobians_up, blitz::Array<double, 4>&)
  %python_attribute(lp_jacobians_db, blitz::Array<double, 4>&)
  %python_attribute(ls_jacobians_db, blitz::Array<double, 3>&)
  %python_attribute(cumtrans, blitz::Array<double, 2>&)
  %python_attribute(lostrans_up, blitz::Array<double, 2>&)
  %python_attribute(lp_cumtrans, blitz::Array<double, 4>&)
  %python_attribute(lp_lostrans_up, blitz::Array<double, 3>&)
  %python_attribute(do_upwelling, bool&)
  %python_attribute(do_dnwelling, bool&)
  %python_attribute(xfine_up, blitz::Array<double, 3>&)
  %python_attribute(wfine_up, blitz::Array<double, 3>&)
  %python_attribute(sunpaths_up, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_up, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_up, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_up, blitz::Array<int, 3>&)
  %python_attribute(xfine_dn, blitz::Array<double, 3>&)
  %python_attribute(wfine_dn, blitz::Array<double, 3>&)
  %python_attribute(sunpaths_dn, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_dn, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_dn, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_dn, blitz::Array<int, 3>&)
  %python_attribute(xfine_up_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_up_p, blitz::Array<double, 3>&)
  %python_attribute(sunpaths_up_p, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_up_p, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_up_p, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_up_p, blitz::Array<int, 3>&)
  %python_attribute(xfine_dn_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_dn_p, blitz::Array<double, 3>&)
  %python_attribute(sunpaths_dn_p, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_dn_p, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_dn_p, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_dn_p, blitz::Array<int, 3>&)
  
  void ss_integral_ilps_dn();
  void ss_integral_ilps_up();
  void ss_integral_ilps_updn();
};


class Fo_Scalarss_Spherfuncs {

public:
  Fo_Scalarss_Spherfuncs(const int& maxmoments_in, const int& maxgeoms_in, const int& nmoments_in, const int& ngeoms_in);
  virtual ~Fo_Scalarss_Spherfuncs();
  std::string print_to_string() const;

  %python_attribute(maxmoments, int&)
  %python_attribute(maxgeoms, int&)
  %python_attribute(nmoments, int&)
  %python_attribute(ngeoms, int&)
  %python_attribute(starter, bool&)
  %python_attribute(do_spherfunc, bool&)
  %python_attribute(cosscat, blitz::Array<double, 1>&)
  %python_attribute(df1, blitz::Array<double, 1>&)
  %python_attribute(df2, blitz::Array<double, 1>&)
  %python_attribute(ss_pleg, blitz::Array<double, 2>&)
  
  void run();
};


class Fo_Thermal_Rtcalcs_I {

public:
  Fo_Thermal_Rtcalcs_I(const int& maxgeoms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& max_user_levels_in, const int& ngeoms_in, const int& nlayers_in, const int& n_user_levels_in, const int& npartials_in);
  virtual ~Fo_Thermal_Rtcalcs_I();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxpartials, int&)
  %python_attribute(maxfine, int&)
  %python_attribute(max_user_levels, int&)
  %python_attribute(do_thermset, bool&)
  %python_attribute(do_deltam_scaling, bool&)
  %python_attribute(do_partials, bool&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(do_sources_dn, blitz::Array<bool, 2>&)
  %python_attribute(do_sources_dn_p, blitz::Array<bool, 2>&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(n_user_levels, int&)
  %python_attribute(user_levels, blitz::Array<int, 1>&)
  %python_attribute(npartials, int&)
  %python_attribute(nfinedivs_p, blitz::Array<int, 2>&)
  %python_attribute(partial_outindex, blitz::Array<int, 1>&)
  %python_attribute(partial_outflag, blitz::Array<bool, 1>&)
  %python_attribute(partial_layeridx, blitz::Array<int, 1>&)
  %python_attribute(bb_input, blitz::Array<double, 1>&)
  %python_attribute(extinction, blitz::Array<double, 1>&)
  %python_attribute(deltaus, blitz::Array<double, 1>&)
  %python_attribute(omega, blitz::Array<double, 1>&)
  %python_attribute(truncfac, blitz::Array<double, 1>&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(losw_paths, blitz::Array<double, 2>&)
  %python_attribute(losp_paths, blitz::Array<double, 2>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(hfine, blitz::Array<double, 3>&)
  %python_attribute(xfine_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_p, blitz::Array<double, 3>&)
  %python_attribute(hfine_p, blitz::Array<double, 3>&)
  %python_attribute(intensity_dta_dn, blitz::Array<double, 2>&)
  %python_attribute(cumsource_dn, blitz::Array<double, 2>&)
  %python_attribute(tcom1, blitz::Array<double, 2>&)
  %python_attribute(do_sources_up, blitz::Array<bool, 2>&)
  %python_attribute(do_sources_up_p, blitz::Array<bool, 2>&)
  %python_attribute(surfbb, double&)
  %python_attribute(user_emissivity, blitz::Array<double, 1>&)
  %python_attribute(intensity_dta_up, blitz::Array<double, 2>&)
  %python_attribute(intensity_dts, blitz::Array<double, 2>&)
  %python_attribute(cumsource_up, blitz::Array<double, 2>&)
  %python_attribute(lostrans_up, blitz::Array<double, 2>&)
  %python_attribute(lostrans_up_p, blitz::Array<double, 2>&)
  %python_attribute(do_upwelling, bool&)
  %python_attribute(do_dnwelling, bool&)
  %python_attribute(xfine_up, blitz::Array<double, 3>&)
  %python_attribute(wfine_up, blitz::Array<double, 3>&)
  %python_attribute(hfine_up, blitz::Array<double, 3>&)
  %python_attribute(xfine_dn, blitz::Array<double, 3>&)
  %python_attribute(wfine_dn, blitz::Array<double, 3>&)
  %python_attribute(hfine_dn, blitz::Array<double, 3>&)
  %python_attribute(xfine_up_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_up_p, blitz::Array<double, 3>&)
  %python_attribute(hfine_up_p, blitz::Array<double, 3>&)
  %python_attribute(xfine_dn_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_dn_p, blitz::Array<double, 3>&)
  %python_attribute(hfine_dn_p, blitz::Array<double, 3>&)
  
  void dte_integral_i_dn();
  void dte_integral_i_up();
  void dte_integral_i_updn();
};


class Fo_Thermal_Rtcalcs_Ilps {

public:
  Fo_Thermal_Rtcalcs_Ilps(const int& maxgeoms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& max_user_levels_in, const int& max_atmoswfs_in, const int& ngeoms_in, const int& nlayers_in, const int& n_user_levels_in, const int& npartials_in, const int& max_surfacewfs_in, const int& n_surfacewfs_in);
  virtual ~Fo_Thermal_Rtcalcs_Ilps();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxpartials, int&)
  %python_attribute(maxfine, int&)
  %python_attribute(max_user_levels, int&)
  %python_attribute(max_atmoswfs, int&)
  %python_attribute(do_thermset, bool&)
  %python_attribute(do_deltam_scaling, bool&)
  %python_attribute(do_partials, bool&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(do_sources_dn, blitz::Array<bool, 2>&)
  %python_attribute(do_sources_dn_p, blitz::Array<bool, 2>&)
  %python_attribute(do_profilewfs, bool&)
  %python_attribute(lvaryflags, blitz::Array<bool, 1>&)
  %python_attribute(lvarynums, blitz::Array<int, 1>&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(n_user_levels, int&)
  %python_attribute(user_levels, blitz::Array<int, 1>&)
  %python_attribute(npartials, int&)
  %python_attribute(nfinedivs_p, blitz::Array<int, 2>&)
  %python_attribute(partial_outindex, blitz::Array<int, 1>&)
  %python_attribute(partial_outflag, blitz::Array<bool, 1>&)
  %python_attribute(partial_layeridx, blitz::Array<int, 1>&)
  %python_attribute(bb_input, blitz::Array<double, 1>&)
  %python_attribute(extinction, blitz::Array<double, 1>&)
  %python_attribute(deltaus, blitz::Array<double, 1>&)
  %python_attribute(omega, blitz::Array<double, 1>&)
  %python_attribute(truncfac, blitz::Array<double, 1>&)
  %python_attribute(l_extinction, blitz::Array<double, 2>&)
  %python_attribute(l_deltaus, blitz::Array<double, 2>&)
  %python_attribute(l_omega, blitz::Array<double, 2>&)
  %python_attribute(l_truncfac, blitz::Array<double, 2>&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(losw_paths, blitz::Array<double, 2>&)
  %python_attribute(losp_paths, blitz::Array<double, 2>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(hfine, blitz::Array<double, 3>&)
  %python_attribute(xfine_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_p, blitz::Array<double, 3>&)
  %python_attribute(hfine_p, blitz::Array<double, 3>&)
  %python_attribute(intensity_dta_dn, blitz::Array<double, 2>&)
  %python_attribute(lp_jacobians_dta_dn, blitz::Array<double, 4>&)
  %python_attribute(tcom1, blitz::Array<double, 2>&)
  %python_attribute(l_tcom1, blitz::Array<double, 3>&)
  %python_attribute(max_surfacewfs, int&)
  %python_attribute(do_sources_up, blitz::Array<bool, 2>&)
  %python_attribute(do_sources_up_p, blitz::Array<bool, 2>&)
  %python_attribute(do_surfacewfs, bool&)
  %python_attribute(n_surfacewfs, int&)
  %python_attribute(surfbb, double&)
  %python_attribute(user_emissivity, blitz::Array<double, 1>&)
  %python_attribute(ls_user_emissivity, blitz::Array<double, 2>&)
  %python_attribute(intensity_dta_up, blitz::Array<double, 2>&)
  %python_attribute(intensity_dts, blitz::Array<double, 2>&)
  %python_attribute(lp_jacobians_dta_up, blitz::Array<double, 4>&)
  %python_attribute(lp_jacobians_dts_up, blitz::Array<double, 4>&)
  %python_attribute(ls_jacobians_dts, blitz::Array<double, 3>&)
  %python_attribute(lostrans_up, blitz::Array<double, 2>&)
  %python_attribute(lostrans_up_p, blitz::Array<double, 2>&)
  %python_attribute(l_lostrans_up, blitz::Array<double, 3>&)
  %python_attribute(l_lostrans_up_p, blitz::Array<double, 3>&)
  %python_attribute(do_upwelling, bool&)
  %python_attribute(do_dnwelling, bool&)
  %python_attribute(xfine_up, blitz::Array<double, 3>&)
  %python_attribute(wfine_up, blitz::Array<double, 3>&)
  %python_attribute(hfine_up, blitz::Array<double, 3>&)
  %python_attribute(xfine_dn, blitz::Array<double, 3>&)
  %python_attribute(wfine_dn, blitz::Array<double, 3>&)
  %python_attribute(hfine_dn, blitz::Array<double, 3>&)
  %python_attribute(xfine_up_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_up_p, blitz::Array<double, 3>&)
  %python_attribute(hfine_up_p, blitz::Array<double, 3>&)
  %python_attribute(xfine_dn_p, blitz::Array<double, 3>&)
  %python_attribute(wfine_dn_p, blitz::Array<double, 3>&)
  %python_attribute(hfine_dn_p, blitz::Array<double, 3>&)
  
  void dte_integral_ilps_dn();
  void dte_integral_ilps_up();
  void dte_integral_ilps_updn();
};

}