// This file was auto-generated

%include "common.i"

%{
#include "first_order_interface.h"
%}



namespace FullPhysics {



class Fo_Dtgeometry_Master {

public:
  Fo_Dtgeometry_Master(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfine_in, const int& ngeoms_in, const int& nlayers_in, const int& nfine_in);
  
  %python_attribute(maxgeoms, int&);
  %python_attribute(maxlayers, int&);
  %python_attribute(maxfine, int&);
  %python_attribute(do_planpar, bool&);
  %python_attribute(do_enhanced_ps, bool&);
  %python_attribute(ngeoms, int&);
  %python_attribute(nlayers, int&);
  %python_attribute(nfine, int&);
  %python_attribute(dtr, double&);
  %python_attribute(eradius, double&);
  %python_attribute(heights, blitz::Array<double, 1>&);
  %python_attribute(alpha_boa, blitz::Array<double, 1>&);
  %python_attribute(donadir, blitz::Array<bool, 1>&);
  %python_attribute(docrit, bool&);
  %python_attribute(acrit, double&);
  %python_attribute(extinc, blitz::Array<double, 1>&);
  %python_attribute(raycon, blitz::Array<double, 1>&);
  %python_attribute(radii, blitz::Array<double, 1>&);
  %python_attribute(alpha, blitz::Array<double, 2>&);
  %python_attribute(cota, blitz::Array<double, 2>&);
  %python_attribute(nfinedivs, blitz::Array<int, 2>&);
  %python_attribute(xfine, blitz::Array<double, 3>&);
  %python_attribute(wfine, blitz::Array<double, 3>&);
  %python_attribute(csqfine, blitz::Array<double, 3>&);
  %python_attribute(cotfine, blitz::Array<double, 3>&);
  %python_attribute(alphafine, blitz::Array<double, 3>&);
  %python_attribute(radiifine, blitz::Array<double, 3>&);
  %python_attribute(ncrit, blitz::Array<int, 1>&);
  %python_attribute(radcrit, blitz::Array<double, 1>&);
  %python_attribute(cotcrit, blitz::Array<double, 1>&);
  %python_attribute(mu1, blitz::Array<double, 1>&);
  %python_attribute(fail, bool&);
  std::string message() const;
  std::string trace() const;
  
  void run();
};


class Fo_Ssgeometry_Master {

public:
  Fo_Ssgeometry_Master(const int& maxgeoms_in, const int& maxszas_in, const int& maxvzas_in, const int& maxazms_in, const int& maxlayers_in, const int& maxfine_in, const int& ngeoms_in, const int& nszas_in, const int& nvzas_in, const int& nazms_in, const int& nlayers_in, const int& nfine_in);
  
  %python_attribute(maxgeoms, int&);
  %python_attribute(maxszas, int&);
  %python_attribute(maxvzas, int&);
  %python_attribute(maxazms, int&);
  %python_attribute(maxlayers, int&);
  %python_attribute(maxfine, int&);
  %python_attribute(do_obsgeom, bool&);
  %python_attribute(do_chapman, bool&);
  %python_attribute(do_planpar, bool&);
  %python_attribute(do_enhanced_ps, bool&);
  %python_attribute(ngeoms, int&);
  %python_attribute(nszas, int&);
  %python_attribute(nvzas, int&);
  %python_attribute(nazms, int&);
  %python_attribute(nlayers, int&);
  %python_attribute(nfine, int&);
  %python_attribute(dtr, double&);
  %python_attribute(pie, double&);
  %python_attribute(vsign, double&);
  %python_attribute(eradius, double&);
  %python_attribute(heights, blitz::Array<double, 1>&);
  %python_attribute(obsgeom_boa, blitz::Array<double, 2>&);
  %python_attribute(alpha_boa, blitz::Array<double, 1>&);
  %python_attribute(theta_boa, blitz::Array<double, 1>&);
  %python_attribute(phi_boa, blitz::Array<double, 1>&);
  %python_attribute(donadir, blitz::Array<bool, 1>&);
  %python_attribute(docrit, bool&);
  %python_attribute(acrit, double&);
  %python_attribute(extinc, blitz::Array<double, 1>&);
  %python_attribute(raycon, blitz::Array<double, 1>&);
  %python_attribute(radii, blitz::Array<double, 1>&);
  %python_attribute(alpha, blitz::Array<double, 2>&);
  %python_attribute(cota, blitz::Array<double, 2>&);
  %python_attribute(nfinedivs, blitz::Array<int, 2>&);
  %python_attribute(xfine, blitz::Array<double, 3>&);
  %python_attribute(wfine, blitz::Array<double, 3>&);
  %python_attribute(csqfine, blitz::Array<double, 3>&);
  %python_attribute(cotfine, blitz::Array<double, 3>&);
  %python_attribute(alphafine, blitz::Array<double, 3>&);
  %python_attribute(radiifine, blitz::Array<double, 3>&);
  %python_attribute(ncrit, blitz::Array<int, 1>&);
  %python_attribute(radcrit, blitz::Array<double, 1>&);
  %python_attribute(cotcrit, blitz::Array<double, 1>&);
  %python_attribute(mu0, blitz::Array<double, 1>&);
  %python_attribute(mu1, blitz::Array<double, 1>&);
  %python_attribute(cosscat, blitz::Array<double, 1>&);
  %python_attribute(chapfacs, blitz::Array<double, 3>&);
  %python_attribute(sunpaths, blitz::Array<double, 3>&);
  %python_attribute(ntraverse, blitz::Array<int, 2>&);
  %python_attribute(sunpathsfine, blitz::Array<double, 4>&);
  %python_attribute(ntraversefine, blitz::Array<int, 3>&);
  %python_attribute(fail, bool&);
  std::string message() const;
  std::string trace() const;
  
  void run();
};

}