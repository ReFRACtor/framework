// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
// This file was auto-generated

// This file was auto-generated

%include "fp_common.i"

%{
#include "first_order_interface.h"
%}

  // MANUAL CHANGE
%base_import(generic_object)
  // MANUAL CHANGE


%fp_shared_ptr(FullPhysics::Fo_Dtgeometry_Master);
%fp_shared_ptr(FullPhysics::Fo_Ssgeometry_Master);
%fp_shared_ptr(FullPhysics::Fo_Scalarss_Rtcalcs_I);
%fp_shared_ptr(FullPhysics::Fo_Scalarss_Rtcalcs_Ilps);
%fp_shared_ptr(FullPhysics::Fo_Scalarss_Spherfuncs);
%fp_shared_ptr(FullPhysics::Fo_Thermal_Rtcalcs_I);
%fp_shared_ptr(FullPhysics::Fo_Thermal_Rtcalcs_Ilpsb);

namespace FullPhysics {



class Fo_Dtgeometry_Master {

public:
  Fo_Dtgeometry_Master(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfine_in, const int& ngeoms_in, const int& nlayers_in, const int& nfine_in);
  virtual ~Fo_Dtgeometry_Master();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxfine, int&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(nfine, int&)
  %python_attribute(dtr, double&)
  %python_attribute(eradius, double&)
  %python_attribute(heights, blitz::Array<double, 1>&)
  %python_attribute(alpha_boa, blitz::Array<double, 1>&)
  %python_attribute(donadir, blitz::Array<bool, 1>&)
  %python_attribute(docrit, bool&)
  %python_attribute(acrit, double&)
  %python_attribute(extinc, blitz::Array<double, 1>&)
  %python_attribute(raycon, blitz::Array<double, 1>&)
  %python_attribute(radii, blitz::Array<double, 1>&)
  %python_attribute(alpha, blitz::Array<double, 2>&)
  %python_attribute(cota, blitz::Array<double, 2>&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(csqfine, blitz::Array<double, 3>&)
  %python_attribute(cotfine, blitz::Array<double, 3>&)
  %python_attribute(alphafine, blitz::Array<double, 3>&)
  %python_attribute(radiifine, blitz::Array<double, 3>&)
  %python_attribute(ncrit, blitz::Array<int, 1>&)
  %python_attribute(radcrit, blitz::Array<double, 1>&)
  %python_attribute(cotcrit, blitz::Array<double, 1>&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(fail, bool&)
  std::string message() const;
  std::string trace() const;
  
  void run();
};


  // MANUAL CHANGE
class Fo_Ssgeometry_Master : public GenericObject {
  // MANUAL CHANGE

public:
  Fo_Ssgeometry_Master(const int& maxgeoms_in, const int& maxszas_in, const int& maxvzas_in, const int& maxazms_in, const int& maxlayers_in, const int& maxfine_in, const int& ngeoms_in, const int& nszas_in, const int& nvzas_in, const int& nazms_in, const int& nlayers_in, const int& nfine_in);
  virtual ~Fo_Ssgeometry_Master();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxszas, int&)
  %python_attribute(maxvzas, int&)
  %python_attribute(maxazms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxfine, int&)
  %python_attribute(do_obsgeom, bool&)
  %python_attribute(do_chapman, bool&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nszas, int&)
  %python_attribute(nvzas, int&)
  %python_attribute(nazms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(nfine, int&)
  %python_attribute(dtr, double&)
  %python_attribute(pie, double&)
  %python_attribute(vsign, double&)
  %python_attribute(eradius, double&)
  %python_attribute(heights, blitz::Array<double, 1>&)
  %python_attribute(obsgeom_boa, blitz::Array<double, 2>&)
  %python_attribute(alpha_boa, blitz::Array<double, 1>&)
  %python_attribute(theta_boa, blitz::Array<double, 1>&)
  %python_attribute(phi_boa, blitz::Array<double, 1>&)
  %python_attribute(donadir, blitz::Array<bool, 1>&)
  %python_attribute(docrit, bool&)
  %python_attribute(acrit, double&)
  %python_attribute(extinc, blitz::Array<double, 1>&)
  %python_attribute(raycon, blitz::Array<double, 1>&)
  %python_attribute(radii, blitz::Array<double, 1>&)
  %python_attribute(alpha, blitz::Array<double, 2>&)
  %python_attribute(cota, blitz::Array<double, 2>&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(csqfine, blitz::Array<double, 3>&)
  %python_attribute(cotfine, blitz::Array<double, 3>&)
  %python_attribute(alphafine, blitz::Array<double, 3>&)
  %python_attribute(radiifine, blitz::Array<double, 3>&)
  %python_attribute(ncrit, blitz::Array<int, 1>&)
  %python_attribute(radcrit, blitz::Array<double, 1>&)
  %python_attribute(cotcrit, blitz::Array<double, 1>&)
  %python_attribute(mu0, blitz::Array<double, 1>&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(cosscat, blitz::Array<double, 1>&)
  %python_attribute(chapfacs, blitz::Array<double, 3>&)
  %python_attribute(sunpaths, blitz::Array<double, 3>&)
  %python_attribute(ntraverse, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine, blitz::Array<int, 3>&)
  %python_attribute(fail, bool&)
  std::string message() const;
  std::string trace() const;
  
  void run();
  // MANUAL CHANGE
  %pickle_serialization();
  // MANUAL CHANGE
};


class Fo_Scalarss_Rtcalcs_I {

public:
  Fo_Scalarss_Rtcalcs_I(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfine_in, const int& ngeoms_in, const int& nlayers_in);
  virtual ~Fo_Scalarss_Rtcalcs_I();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxfine, int&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_regular_ps, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(donadir, blitz::Array<bool, 1>&)
  %python_attribute(do_sleave, bool&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(aclevel, int&)
  %python_attribute(reflec, blitz::Array<double, 1>&)
  %python_attribute(slterm, blitz::Array<double, 1>&)
  %python_attribute(extinction, blitz::Array<double, 1>&)
  %python_attribute(deltaus, blitz::Array<double, 1>&)
  %python_attribute(exactscat_up, blitz::Array<double, 2>&)
  %python_attribute(flux, double&)
  %python_attribute(mu0, blitz::Array<double, 1>&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(ncrit, blitz::Array<int, 1>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(csqfine, blitz::Array<double, 3>&)
  %python_attribute(cotfine, blitz::Array<double, 3>&)
  %python_attribute(raycon, blitz::Array<double, 1>&)
  %python_attribute(cota, blitz::Array<double, 2>&)
  %python_attribute(sunpaths, blitz::Array<double, 3>&)
  %python_attribute(ntraverse, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine, blitz::Array<int, 3>&)
  %python_attribute(intensity_up, blitz::Array<double, 1>&)
  %python_attribute(intensity_db, blitz::Array<double, 1>&)
  %python_attribute(cumsource_up, blitz::Array<double, 2>&)
  
  void ss_integral_i_up();
};


  // MANUAL CHANGE
class Fo_Scalarss_Rtcalcs_Ilps : public GenericObject {
  // MANUAL CHANGE

public:
  Fo_Scalarss_Rtcalcs_Ilps(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfine_in, const int& max_atmoswfs_in, const int& max_surfacewfs_in, const int& ngeoms_in, const int& nlayers_in);
  virtual ~Fo_Scalarss_Rtcalcs_Ilps();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxfine, int&)
  %python_attribute(max_atmoswfs, int&)
  %python_attribute(max_surfacewfs, int&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_regular_ps, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(donadir, blitz::Array<bool, 1>&)
  %python_attribute(do_sleave, bool&)
  %python_attribute(do_profilewfs, bool&)
  %python_attribute(do_reflecwfs, bool&)
  %python_attribute(do_sleavewfs, bool&)
  %python_attribute(lvaryflags, blitz::Array<bool, 1>&)
  %python_attribute(lvarynums, blitz::Array<int, 1>&)
  %python_attribute(n_reflecwfs, int&)
  %python_attribute(n_sleavewfs, int&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(aclevel, int&)
  %python_attribute(reflec, blitz::Array<double, 1>&)
  %python_attribute(slterm, blitz::Array<double, 1>&)
  %python_attribute(extinction, blitz::Array<double, 1>&)
  %python_attribute(deltaus, blitz::Array<double, 1>&)
  %python_attribute(exactscat_up, blitz::Array<double, 2>&)
  %python_attribute(flux, double&)
  %python_attribute(ls_reflec, blitz::Array<double, 2>&)
  %python_attribute(lssl_slterm, blitz::Array<double, 2>&)
  %python_attribute(l_extinction, blitz::Array<double, 2>&)
  %python_attribute(l_deltaus, blitz::Array<double, 2>&)
  %python_attribute(l_exactscat_up, blitz::Array<double, 3>&)
  %python_attribute(mu0, blitz::Array<double, 1>&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(ncrit, blitz::Array<int, 1>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(csqfine, blitz::Array<double, 3>&)
  %python_attribute(cotfine, blitz::Array<double, 3>&)
  %python_attribute(raycon, blitz::Array<double, 1>&)
  %python_attribute(cota, blitz::Array<double, 2>&)
  %python_attribute(sunpaths, blitz::Array<double, 3>&)
  %python_attribute(ntraverse, blitz::Array<int, 2>&)
  %python_attribute(sunpaths_fine, blitz::Array<double, 4>&)
  %python_attribute(ntraverse_fine, blitz::Array<int, 3>&)
  %python_attribute(intensity_up, blitz::Array<double, 1>&)
  %python_attribute(intensity_db, blitz::Array<double, 1>&)
  %python_attribute(lp_jacobians_up, blitz::Array<double, 3>&)
  %python_attribute(lp_jacobians_db, blitz::Array<double, 3>&)
  %python_attribute(ls_jacobians_db, blitz::Array<double, 2>&)
  
  void ss_integral_ilps_up();
  // MANUAL CHANGE
  %pickle_serialization();
  // MANUAL CHANGE
};


  // MANUAL CHANGE
class Fo_Scalarss_Spherfuncs : public GenericObject {
  // MANUAL CHANGE

public:
  Fo_Scalarss_Spherfuncs(const bool& starter_in, const int& maxmoms_in, const int& maxgeoms_in, const int& nmoms_in, const int& ngeoms_in);
  virtual ~Fo_Scalarss_Spherfuncs();
  std::string print_to_string() const;

  %python_attribute(starter, bool&)
  %python_attribute(maxmoms, int&)
  %python_attribute(maxgeoms, int&)
  %python_attribute(nmoms, int&)
  %python_attribute(ngeoms, int&)
  %python_attribute(df1, blitz::Array<double, 1>&)
  %python_attribute(df2, blitz::Array<double, 1>&)
  %python_attribute(cosscat, blitz::Array<double, 1>&)
  %python_attribute(ss_pleg, blitz::Array<double, 2>&)
  
  void run();
  // MANUAL CHANGE
  %pickle_serialization();
  // MANUAL CHANGE
};


class Fo_Thermal_Rtcalcs_I {

public:
  Fo_Thermal_Rtcalcs_I(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfinelayers_in, const int& max_user_levels_in, const int& ngeoms_in, const int& nlayers_in, const int& n_user_levels_in);
  virtual ~Fo_Thermal_Rtcalcs_I();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxfinelayers, int&)
  %python_attribute(max_user_levels, int&)
  %python_attribute(do_thermset, bool&)
  %python_attribute(do_deltam_scaling, bool&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_regular_ps, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(donadir, blitz::Array<bool, 1>&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(n_user_levels, int&)
  %python_attribute(user_levels, blitz::Array<int, 1>&)
  %python_attribute(bb_input, blitz::Array<double, 1>&)
  %python_attribute(extinction, blitz::Array<double, 1>&)
  %python_attribute(deltaus, blitz::Array<double, 1>&)
  %python_attribute(omega, blitz::Array<double, 1>&)
  %python_attribute(truncfac, blitz::Array<double, 1>&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(ncrit, blitz::Array<int, 1>&)
  %python_attribute(radcrit, blitz::Array<double, 1>&)
  %python_attribute(cotcrit, blitz::Array<double, 1>&)
  %python_attribute(raycon, blitz::Array<double, 1>&)
  %python_attribute(radii, blitz::Array<double, 1>&)
  %python_attribute(cota, blitz::Array<double, 2>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(csqfine, blitz::Array<double, 3>&)
  %python_attribute(cotfine, blitz::Array<double, 3>&)
  %python_attribute(intensity_dta_dn, blitz::Array<double, 2>&)
  %python_attribute(cumsource_dn, blitz::Array<double, 2>&)
  %python_attribute(tcom1, blitz::Array<double, 2>&)
  %python_attribute(surfbb, double&)
  %python_attribute(user_emissivity, blitz::Array<double, 1>&)
  %python_attribute(intensity_dta_up, blitz::Array<double, 2>&)
  %python_attribute(intensity_dts, blitz::Array<double, 2>&)
  %python_attribute(cumsource_up, blitz::Array<double, 2>&)
  %python_attribute(do_upwelling, bool&)
  %python_attribute(do_dnwelling, bool&)
  
  void dte_integral_i_dn();
  void dte_integral_i_up();
  void dte_integral_i_updn();
};


class Fo_Thermal_Rtcalcs_Ilpsb {

public:
  Fo_Thermal_Rtcalcs_Ilpsb(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfinelayers_in, const int& max_user_levels_in, const int& max_atmoswfs_in, const int& ngeoms_in, const int& nlayers_in, const int& n_user_levels_in, const int& max_surfacewfs_in);
  virtual ~Fo_Thermal_Rtcalcs_Ilpsb();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxfinelayers, int&)
  %python_attribute(max_user_levels, int&)
  %python_attribute(max_atmoswfs, int&)
  %python_attribute(do_abbwf, bool&)
  %python_attribute(do_thermset, bool&)
  %python_attribute(do_deltam_scaling, bool&)
  %python_attribute(do_planpar, bool&)
  %python_attribute(do_regular_ps, bool&)
  %python_attribute(do_enhanced_ps, bool&)
  %python_attribute(donadir, blitz::Array<bool, 1>&)
  %python_attribute(do_profilewfs, bool&)
  %python_attribute(lvaryflags, blitz::Array<bool, 1>&)
  %python_attribute(lvarynums, blitz::Array<int, 1>&)
  %python_attribute(ngeoms, int&)
  %python_attribute(nlayers, int&)
  %python_attribute(nfinedivs, blitz::Array<int, 2>&)
  %python_attribute(n_user_levels, int&)
  %python_attribute(user_levels, blitz::Array<int, 1>&)
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
  %python_attribute(ncrit, blitz::Array<int, 1>&)
  %python_attribute(radcrit, blitz::Array<double, 1>&)
  %python_attribute(cotcrit, blitz::Array<double, 1>&)
  %python_attribute(raycon, blitz::Array<double, 1>&)
  %python_attribute(cota, blitz::Array<double, 2>&)
  %python_attribute(radii, blitz::Array<double, 1>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(csqfine, blitz::Array<double, 3>&)
  %python_attribute(cotfine, blitz::Array<double, 3>&)
  %python_attribute(intensity_dta_dn, blitz::Array<double, 2>&)
  %python_attribute(lp_jacobians_dta_dn, blitz::Array<double, 4>&)
  %python_attribute(lab_jacobians_dta_dn, blitz::Array<double, 3>&)
  %python_attribute(tcom, blitz::Array<double, 2>&)
  %python_attribute(l_tcom, blitz::Array<double, 3>&)
  %python_attribute(lb_tcom1, blitz::Array<double, 2>&)
  %python_attribute(lb_tcom2, blitz::Array<double, 2>&)
  %python_attribute(max_surfacewfs, int&)
  %python_attribute(do_sbbwf, bool&)
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
  %python_attribute(lab_jacobians_dta_up, blitz::Array<double, 3>&)
  %python_attribute(lsb_jacobians_dts, blitz::Array<double, 2>&)
  %python_attribute(do_upwelling, bool&)
  %python_attribute(do_dnwelling, bool&)
  
  void dte_integral_ilpsb_dn();
  void dte_integral_ilpsb_up();
  void dte_integral_ilpsb_updn();
};

}
