// This file was auto-generated

%include "common.i"

%{
#include "first_order_interface.h"
%}



%fp_shared_ptr(FullPhysics::Fo_Dtgeometry_Master);
%fp_shared_ptr(FullPhysics::Fo_Ssgeometry_Master);
%fp_shared_ptr(FullPhysics::Fo_Scalarss_Rtcalcs_I);
%fp_shared_ptr(FullPhysics::Fo_Scalarss_Spherfuncs);

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


class Fo_Ssgeometry_Master {

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
};


class Fo_Scalarss_Rtcalcs_I {

public:
  Fo_Scalarss_Rtcalcs_I(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfine_in, const int& max_user_levels_in, const int& ngeoms_in, const int& nlayers_in, const int& n_user_levels_in);
  virtual ~Fo_Scalarss_Rtcalcs_I();
  std::string print_to_string() const;

  %python_attribute(maxgeoms, int&)
  %python_attribute(maxlayers, int&)
  %python_attribute(maxfine, int&)
  %python_attribute(max_user_levels, int&)
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
  %python_attribute(extinction, blitz::Array<double, 1>&)
  %python_attribute(deltaus, blitz::Array<double, 1>&)
  %python_attribute(exactscat_dn, blitz::Array<double, 2>&)
  %python_attribute(flux, double&)
  %python_attribute(mu1, blitz::Array<double, 1>&)
  %python_attribute(ncrit, blitz::Array<int, 1>&)
  %python_attribute(radcrit, blitz::Array<double, 1>&)
  %python_attribute(cotcrit, blitz::Array<double, 1>&)
  %python_attribute(xfine, blitz::Array<double, 3>&)
  %python_attribute(wfine, blitz::Array<double, 3>&)
  %python_attribute(csqfine, blitz::Array<double, 3>&)
  %python_attribute(cotfine, blitz::Array<double, 3>&)
  %python_attribute(raycon, blitz::Array<double, 1>&)
  %python_attribute(radii, blitz::Array<double, 1>&)
  %python_attribute(cota, blitz::Array<double, 2>&)
  %python_attribute(sunpaths, blitz::Array<double, 3>&)
  %python_attribute(ntraverse, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine, blitz::Array<int, 3>&)
  %python_attribute(intensity_dn, blitz::Array<double, 2>&)
  %python_attribute(cumsource_dn, blitz::Array<double, 2>&)
  %python_attribute(reflec, blitz::Array<double, 1>&)
  %python_attribute(exactscat_up, blitz::Array<double, 2>&)
  %python_attribute(mu0, blitz::Array<double, 1>&)
  %python_attribute(intensity_up, blitz::Array<double, 2>&)
  %python_attribute(intensity_db, blitz::Array<double, 2>&)
  %python_attribute(cumsource_up, blitz::Array<double, 2>&)
  %python_attribute(do_upwelling, bool&)
  %python_attribute(do_dnwelling, bool&)
  %python_attribute(sunpaths_up, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_up, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_up, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_up, blitz::Array<int, 3>&)
  %python_attribute(sunpaths_dn, blitz::Array<double, 3>&)
  %python_attribute(ntraverse_dn, blitz::Array<int, 2>&)
  %python_attribute(sunpathsfine_dn, blitz::Array<double, 4>&)
  %python_attribute(ntraversefine_dn, blitz::Array<int, 3>&)
  
  void ss_integral_i_dn();
  void ss_integral_i_up();
  void ss_integral_i_updn();
};


class Fo_Scalarss_Spherfuncs {

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
};

}