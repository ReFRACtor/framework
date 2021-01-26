#ifndef FIRST_ORDER_INTERFACE_H
#define FIRST_ORDER_INTERFACE_H

#include <iostream>
#include <blitz/array.h>

#include "fp_exception.h"

#define FIRST_ORDER_DEBUG 0

/* This file was auto-generated */

namespace FullPhysics {

//-----------------------------------------------------------------------
// Links to module: "fo_dtgeometry_master_m" in file: "FO_DTgeometry_master.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_dtgeometry_master_m_fo_dtgeometry_master_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxfine_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const int* ngeoms_in, const int* nlayers_in, const int* nfine_in, const double* dtr_in, const double* eradius_in, const double* heights_in, const double* alpha_boa_in, const bool* donadir_in, const bool* docrit_in, const double* acrit_in, const double* extinc_in, const double* raycon_in, const double* radii_in, const double* alpha_in, const double* cota_in, const int* nfinedivs_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* alphafine_in, const double* radiifine_in, const int* ncrit_in, const double* radcrit_in, const double* cotcrit_in, const double* mu1_in, const bool* fail_in, const int* message_len, const char* message_in, const int* trace_len, const char* trace_in);
}

class Fo_Dtgeometry_Master : public virtual GenericObject {

public:
  Fo_Dtgeometry_Master(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfine_in, const int& ngeoms_in, const int& nlayers_in, const int& nfine_in) : maxgeoms_(maxgeoms_in), maxlayers_(maxlayers_in), maxfine_(maxfine_in), ngeoms_(ngeoms_in), nlayers_(nlayers_in), nfine_(nfine_in) 
  { 
    do_planpar_ = false;
    do_enhanced_ps_ = false;
    dtr_ = 0;
    eradius_ = 0;
    heights_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    heights_ = 0;
    alpha_boa_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    alpha_boa_ = 0;
    donadir_.reference( blitz::Array<bool, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    donadir_ = false;
    docrit_ = false;
    acrit_ = 0;
    extinc_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    extinc_ = 0;
    raycon_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    raycon_ = 0;
    radii_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    radii_ = 0;
    alpha_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    alpha_ = 0;
    cota_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cota_ = 0;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    csqfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    csqfine_ = 0;
    cotfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    cotfine_ = 0;
    alphafine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    alphafine_ = 0;
    radiifine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    radiifine_ = 0;
    ncrit_.reference( blitz::Array<int, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    ncrit_ = 0;
    radcrit_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    radcrit_ = 0;
    cotcrit_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    cotcrit_ = 0;
    mu1_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu1_ = 0;
    fail_ = false;
    message_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    message_ = '\0';
    trace_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    trace_ = '\0';
    // Initialize type pointers
    
  }

  virtual ~Fo_Dtgeometry_Master() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxfine() const {
    return maxfine_;
  }

  

  const bool& do_planpar() const {
    return do_planpar_;
  }

  void do_planpar(const bool& do_planpar_in) {
    do_planpar_ = do_planpar_in;
  }

  

  const bool& do_enhanced_ps() const {
    return do_enhanced_ps_;
  }

  void do_enhanced_ps(const bool& do_enhanced_ps_in) {
    do_enhanced_ps_ = do_enhanced_ps_in;
  }

  

  const int& ngeoms() const {
    return ngeoms_;
  }

  

  const int& nlayers() const {
    return nlayers_;
  }

  

  const int& nfine() const {
    return nfine_;
  }

  

  const double& dtr() const {
    return dtr_;
  }

  void dtr(const double& dtr_in) {
    dtr_ = dtr_in;
  }

  

  const double& eradius() const {
    return eradius_;
  }

  void eradius(const double& eradius_in) {
    eradius_ = eradius_in;
  }

  

  const blitz::Array<double, 1>& heights() const {
    return heights_;
  }

  void heights(const blitz::Array<double, 1>& heights_in) {
    heights_ = heights_in;
  }

  

  const blitz::Array<double, 1>& alpha_boa() const {
    return alpha_boa_;
  }

  void alpha_boa(const blitz::Array<double, 1>& alpha_boa_in) {
    alpha_boa_ = alpha_boa_in;
  }

  

  const blitz::Array<bool, 1>& donadir() const {
    return donadir_;
  }

  void donadir(const blitz::Array<bool, 1>& donadir_in) {
    donadir_ = donadir_in;
  }

  

  const bool& docrit() const {
    return docrit_;
  }

  void docrit(const bool& docrit_in) {
    docrit_ = docrit_in;
  }

  

  const double& acrit() const {
    return acrit_;
  }

  void acrit(const double& acrit_in) {
    acrit_ = acrit_in;
  }

  

  const blitz::Array<double, 1>& extinc() const {
    return extinc_;
  }

  void extinc(const blitz::Array<double, 1>& extinc_in) {
    extinc_ = extinc_in;
  }

  

  const blitz::Array<double, 1>& raycon() const {
    return raycon_;
  }

  void raycon(const blitz::Array<double, 1>& raycon_in) {
    raycon_ = raycon_in;
  }

  

  const blitz::Array<double, 1>& radii() const {
    return radii_;
  }

  void radii(const blitz::Array<double, 1>& radii_in) {
    radii_ = radii_in;
  }

  

  const blitz::Array<double, 2>& alpha() const {
    return alpha_;
  }

  void alpha(const blitz::Array<double, 2>& alpha_in) {
    alpha_ = alpha_in;
  }

  

  const blitz::Array<double, 2>& cota() const {
    return cota_;
  }

  void cota(const blitz::Array<double, 2>& cota_in) {
    cota_ = cota_in;
  }

  

  const blitz::Array<int, 2>& nfinedivs() const {
    return nfinedivs_;
  }

  void nfinedivs(const blitz::Array<int, 2>& nfinedivs_in) {
    nfinedivs_ = nfinedivs_in;
  }

  

  const blitz::Array<double, 3>& xfine() const {
    return xfine_;
  }

  void xfine(const blitz::Array<double, 3>& xfine_in) {
    xfine_ = xfine_in;
  }

  

  const blitz::Array<double, 3>& wfine() const {
    return wfine_;
  }

  void wfine(const blitz::Array<double, 3>& wfine_in) {
    wfine_ = wfine_in;
  }

  

  const blitz::Array<double, 3>& csqfine() const {
    return csqfine_;
  }

  void csqfine(const blitz::Array<double, 3>& csqfine_in) {
    csqfine_ = csqfine_in;
  }

  

  const blitz::Array<double, 3>& cotfine() const {
    return cotfine_;
  }

  void cotfine(const blitz::Array<double, 3>& cotfine_in) {
    cotfine_ = cotfine_in;
  }

  

  const blitz::Array<double, 3>& alphafine() const {
    return alphafine_;
  }

  void alphafine(const blitz::Array<double, 3>& alphafine_in) {
    alphafine_ = alphafine_in;
  }

  

  const blitz::Array<double, 3>& radiifine() const {
    return radiifine_;
  }

  void radiifine(const blitz::Array<double, 3>& radiifine_in) {
    radiifine_ = radiifine_in;
  }

  

  const blitz::Array<int, 1>& ncrit() const {
    return ncrit_;
  }

  

  const blitz::Array<double, 1>& radcrit() const {
    return radcrit_;
  }

  

  const blitz::Array<double, 1>& cotcrit() const {
    return cotcrit_;
  }

  

  const blitz::Array<double, 1>& mu1() const {
    return mu1_;
  }

  

  const bool& fail() const {
    return fail_;
  }

  

  std::string message() const {
    std::string message_ret;
    message_ret = ( std::string(std::string(message_(blitz::Range::all()).begin(), message_(blitz::Range::all()).end()).c_str()) );
    return message_ret;
  }

  

  std::string trace() const {
    std::string trace_ret;
    trace_ret = ( std::string(std::string(trace_(blitz::Range::all()).begin(), trace_(blitz::Range::all()).end()).c_str()) );
    return trace_ret;
  }

  

  
  void run() {
    int message_len = (int) message_.extent(0) - 1;
    int trace_len = (int) trace_.extent(0) - 1;
    
    fo_dtgeometry_master_m_fo_dtgeometry_master_wrap(&maxgeoms_, &maxlayers_, &maxfine_, &do_planpar_, &do_enhanced_ps_, &ngeoms_, &nlayers_, &nfine_, &dtr_, &eradius_, heights_.dataFirst(), alpha_boa_.dataFirst(), donadir_.dataFirst(), &docrit_, &acrit_, extinc_.dataFirst(), raycon_.dataFirst(), radii_.dataFirst(), alpha_.dataFirst(), cota_.dataFirst(), nfinedivs_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), csqfine_.dataFirst(), cotfine_.dataFirst(), alphafine_.dataFirst(), radiifine_.dataFirst(), ncrit_.dataFirst(), radcrit_.dataFirst(), cotcrit_.dataFirst(), mu1_.dataFirst(), &fail_, &message_len, message_.dataFirst(), &trace_len, trace_.dataFirst());
    
  }

  friend std::ostream& operator<<(std::ostream &output_stream, const Fo_Dtgeometry_Master &obj) {
    output_stream << "Fo_Dtgeometry_Master:" << std::endl
      << "      maxgeoms: " << obj.maxgeoms()  << std::endl
      << "     maxlayers: " << obj.maxlayers()  << std::endl
      << "       maxfine: " << obj.maxfine()  << std::endl
      << "    do_planpar: " << obj.do_planpar()  << std::endl
      << "do_enhanced_ps: " << obj.do_enhanced_ps()  << std::endl
      << "        ngeoms: " << obj.ngeoms()  << std::endl
      << "       nlayers: " << obj.nlayers()  << std::endl
      << "         nfine: " << obj.nfine()  << std::endl
      << "           dtr: " << obj.dtr()  << std::endl
      << "       eradius: " << obj.eradius()  << std::endl
      << "       heights: " << std::endl << obj.heights()  << std::endl
      << "     alpha_boa: " << std::endl << obj.alpha_boa()  << std::endl
      << "       donadir: " << std::endl << obj.donadir()  << std::endl
      << "        docrit: " << obj.docrit()  << std::endl
      << "         acrit: " << obj.acrit()  << std::endl
      << "        extinc: " << std::endl << obj.extinc()  << std::endl
      << "        raycon: " << std::endl << obj.raycon()  << std::endl
      << "         radii: " << std::endl << obj.radii()  << std::endl
      << "         alpha: " << std::endl << obj.alpha()  << std::endl
      << "          cota: " << std::endl << obj.cota()  << std::endl
      << "     nfinedivs: " << std::endl << obj.nfinedivs()  << std::endl
      << "         xfine: " << std::endl << obj.xfine()  << std::endl
      << "         wfine: " << std::endl << obj.wfine()  << std::endl
      << "       csqfine: " << std::endl << obj.csqfine()  << std::endl
      << "       cotfine: " << std::endl << obj.cotfine()  << std::endl
      << "     alphafine: " << std::endl << obj.alphafine()  << std::endl
      << "     radiifine: " << std::endl << obj.radiifine()  << std::endl
      << "         ncrit: " << std::endl << obj.ncrit()  << std::endl
      << "       radcrit: " << std::endl << obj.radcrit()  << std::endl
      << "       cotcrit: " << std::endl << obj.cotcrit()  << std::endl
      << "           mu1: " << std::endl << obj.mu1()  << std::endl
      << "          fail: " << obj.fail()  << std::endl
      << "       message: " << "\"" << obj.message() << "\"" << std::endl
      << "         trace: " << "\"" << obj.trace() << "\"" << std::endl;
    return output_stream;

  }

private:
  int maxgeoms_;
  int maxlayers_;
  int maxfine_;
  bool do_planpar_;
  bool do_enhanced_ps_;
  int ngeoms_;
  int nlayers_;
  int nfine_;
  double dtr_;
  double eradius_;
  blitz::Array<double, 1> heights_;
  blitz::Array<double, 1> alpha_boa_;
  blitz::Array<bool, 1> donadir_;
  bool docrit_;
  double acrit_;
  blitz::Array<double, 1> extinc_;
  blitz::Array<double, 1> raycon_;
  blitz::Array<double, 1> radii_;
  blitz::Array<double, 2> alpha_;
  blitz::Array<double, 2> cota_;
  blitz::Array<int, 2> nfinedivs_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> csqfine_;
  blitz::Array<double, 3> cotfine_;
  blitz::Array<double, 3> alphafine_;
  blitz::Array<double, 3> radiifine_;
  blitz::Array<int, 1> ncrit_;
  blitz::Array<double, 1> radcrit_;
  blitz::Array<double, 1> cotcrit_;
  blitz::Array<double, 1> mu1_;
  bool fail_;
  blitz::Array<char, 1> message_;
  blitz::Array<char, 1> trace_;
};

//-----------------------------------------------------------------------
// Links to module: "fo_ssgeometry_master_m" in file: "FO_SSgeometry_master.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_ssgeometry_master_m_fo_ssgeometry_master_wrap(const int* maxgeoms_in, const int* maxszas_in, const int* maxvzas_in, const int* maxazms_in, const int* maxlayers_in, const int* maxfine_in, const bool* do_obsgeom_in, const bool* do_chapman_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const int* ngeoms_in, const int* nszas_in, const int* nvzas_in, const int* nazms_in, const int* nlayers_in, const int* nfine_in, const double* dtr_in, const double* pie_in, const double* vsign_in, const double* eradius_in, const double* heights_in, const double* obsgeom_boa_in, const double* alpha_boa_in, const double* theta_boa_in, const double* phi_boa_in, const bool* donadir_in, const bool* docrit_in, const double* acrit_in, const double* extinc_in, const double* raycon_in, const double* radii_in, const double* alpha_in, const double* cota_in, const int* nfinedivs_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* alphafine_in, const double* radiifine_in, const int* ncrit_in, const double* radcrit_in, const double* cotcrit_in, const double* mu0_in, const double* mu1_in, const double* cosscat_in, const double* chapfacs_in, const double* sunpaths_in, const int* ntraverse_in, const double* sunpathsfine_in, const int* ntraversefine_in, const bool* fail_in, const int* message_len, const char* message_in, const int* trace_len, const char* trace_in);
}

  // MANUAL CHANGE
class Fo_Ssgeometry_Master : public Printable<Fo_Ssgeometry_Master> {
  // MANUAL CHANGE

public:
  Fo_Ssgeometry_Master(const int& maxgeoms_in, const int& maxszas_in, const int& maxvzas_in, const int& maxazms_in, const int& maxlayers_in, const int& maxfine_in, const int& ngeoms_in, const int& nszas_in, const int& nvzas_in, const int& nazms_in, const int& nlayers_in, const int& nfine_in) : maxgeoms_(maxgeoms_in), maxszas_(maxszas_in), maxvzas_(maxvzas_in), maxazms_(maxazms_in), maxlayers_(maxlayers_in), maxfine_(maxfine_in), ngeoms_(ngeoms_in), nszas_(nszas_in), nvzas_(nvzas_in), nazms_(nazms_in), nlayers_(nlayers_in), nfine_(nfine_in) 
  { 
    do_obsgeom_ = false;
    do_chapman_ = false;
    do_planpar_ = false;
    do_enhanced_ps_ = false;
    dtr_ = 0;
    pie_ = 0;
    vsign_ = 0;
    eradius_ = 0;
    heights_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    heights_ = 0;
    obsgeom_boa_.reference( blitz::Array<double, 2>(maxgeoms_, 3, blitz::ColumnMajorArray<2>()) );
    obsgeom_boa_ = 0;
    alpha_boa_.reference( blitz::Array<double, 1>(maxvzas_, blitz::ColumnMajorArray<1>()) );
    alpha_boa_ = 0;
    theta_boa_.reference( blitz::Array<double, 1>(maxszas_, blitz::ColumnMajorArray<1>()) );
    theta_boa_ = 0;
    phi_boa_.reference( blitz::Array<double, 1>(maxazms_, blitz::ColumnMajorArray<1>()) );
    phi_boa_ = 0;
    donadir_.reference( blitz::Array<bool, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    donadir_ = false;
    docrit_ = false;
    acrit_ = 0;
    extinc_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    extinc_ = 0;
    raycon_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    raycon_ = 0;
    radii_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    radii_ = 0;
    alpha_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    alpha_ = 0;
    cota_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cota_ = 0;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    csqfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    csqfine_ = 0;
    cotfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    cotfine_ = 0;
    alphafine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    alphafine_ = 0;
    radiifine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    radiifine_ = 0;
    ncrit_.reference( blitz::Array<int, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    ncrit_ = 0;
    radcrit_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    radcrit_ = 0;
    cotcrit_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    cotcrit_ = 0;
    mu0_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu0_ = 0;
    mu1_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu1_ = 0;
    cosscat_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    cosscat_ = 0;
    chapfacs_.reference( blitz::Array<double, 3>(maxlayers_, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    chapfacs_ = 0;
    sunpaths_.reference( blitz::Array<double, 3>(maxlayers_-0+1, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_ = 0;
    ntraverse_.reference( blitz::Array<int, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_ = 0;
    sunpathsfine_.reference( blitz::Array<double, 4>(maxlayers_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_ = 0;
    ntraversefine_.reference( blitz::Array<int, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_ = 0;
    fail_ = false;
    message_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    message_ = '\0';
    trace_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    trace_ = '\0';
    // Initialize type pointers
    
  }

  virtual ~Fo_Ssgeometry_Master() = default;

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& maxszas() const {
    return maxszas_;
  }

  

  const int& maxvzas() const {
    return maxvzas_;
  }

  

  const int& maxazms() const {
    return maxazms_;
  }

  

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxfine() const {
    return maxfine_;
  }

  

  const bool& do_obsgeom() const {
    return do_obsgeom_;
  }

  void do_obsgeom(const bool& do_obsgeom_in) {
    do_obsgeom_ = do_obsgeom_in;
  }

  

  const bool& do_chapman() const {
    return do_chapman_;
  }

  void do_chapman(const bool& do_chapman_in) {
    do_chapman_ = do_chapman_in;
  }

  

  const bool& do_planpar() const {
    return do_planpar_;
  }

  void do_planpar(const bool& do_planpar_in) {
    do_planpar_ = do_planpar_in;
  }

  

  const bool& do_enhanced_ps() const {
    return do_enhanced_ps_;
  }

  void do_enhanced_ps(const bool& do_enhanced_ps_in) {
    do_enhanced_ps_ = do_enhanced_ps_in;
  }

  

  const int& ngeoms() const {
    return ngeoms_;
  }

  

  const int& nszas() const {
    return nszas_;
  }

  

  const int& nvzas() const {
    return nvzas_;
  }

  

  const int& nazms() const {
    return nazms_;
  }

  

  const int& nlayers() const {
    return nlayers_;
  }

  

  const int& nfine() const {
    return nfine_;
  }

  

  const double& dtr() const {
    return dtr_;
  }

  void dtr(const double& dtr_in) {
    dtr_ = dtr_in;
  }

  

  const double& pie() const {
    return pie_;
  }

  void pie(const double& pie_in) {
    pie_ = pie_in;
  }

  

  const double& vsign() const {
    return vsign_;
  }

  void vsign(const double& vsign_in) {
    vsign_ = vsign_in;
  }

  

  const double& eradius() const {
    return eradius_;
  }

  void eradius(const double& eradius_in) {
    eradius_ = eradius_in;
  }

  

  const blitz::Array<double, 1>& heights() const {
    return heights_;
  }

  void heights(const blitz::Array<double, 1>& heights_in) {
    heights_ = heights_in;
  }

  

  const blitz::Array<double, 2>& obsgeom_boa() const {
    return obsgeom_boa_;
  }

  void obsgeom_boa(const blitz::Array<double, 2>& obsgeom_boa_in) {
    obsgeom_boa_ = obsgeom_boa_in;
  }

  

  const blitz::Array<double, 1>& alpha_boa() const {
    return alpha_boa_;
  }

  void alpha_boa(const blitz::Array<double, 1>& alpha_boa_in) {
    alpha_boa_ = alpha_boa_in;
  }

  

  const blitz::Array<double, 1>& theta_boa() const {
    return theta_boa_;
  }

  void theta_boa(const blitz::Array<double, 1>& theta_boa_in) {
    theta_boa_ = theta_boa_in;
  }

  

  const blitz::Array<double, 1>& phi_boa() const {
    return phi_boa_;
  }

  void phi_boa(const blitz::Array<double, 1>& phi_boa_in) {
    phi_boa_ = phi_boa_in;
  }

  

  const blitz::Array<bool, 1>& donadir() const {
    return donadir_;
  }

  void donadir(const blitz::Array<bool, 1>& donadir_in) {
    donadir_ = donadir_in;
  }

  

  const bool& docrit() const {
    return docrit_;
  }

  void docrit(const bool& docrit_in) {
    docrit_ = docrit_in;
  }

  

  const double& acrit() const {
    return acrit_;
  }

  void acrit(const double& acrit_in) {
    acrit_ = acrit_in;
  }

  

  const blitz::Array<double, 1>& extinc() const {
    return extinc_;
  }

  void extinc(const blitz::Array<double, 1>& extinc_in) {
    extinc_ = extinc_in;
  }

  

  const blitz::Array<double, 1>& raycon() const {
    return raycon_;
  }

  void raycon(const blitz::Array<double, 1>& raycon_in) {
    raycon_ = raycon_in;
  }

  

  const blitz::Array<double, 1>& radii() const {
    return radii_;
  }

  void radii(const blitz::Array<double, 1>& radii_in) {
    radii_ = radii_in;
  }

  

  const blitz::Array<double, 2>& alpha() const {
    return alpha_;
  }

  void alpha(const blitz::Array<double, 2>& alpha_in) {
    alpha_ = alpha_in;
  }

  

  const blitz::Array<double, 2>& cota() const {
    return cota_;
  }

  void cota(const blitz::Array<double, 2>& cota_in) {
    cota_ = cota_in;
  }

  

  const blitz::Array<int, 2>& nfinedivs() const {
    return nfinedivs_;
  }

  void nfinedivs(const blitz::Array<int, 2>& nfinedivs_in) {
    nfinedivs_ = nfinedivs_in;
  }

  

  const blitz::Array<double, 3>& xfine() const {
    return xfine_;
  }

  void xfine(const blitz::Array<double, 3>& xfine_in) {
    xfine_ = xfine_in;
  }

  

  const blitz::Array<double, 3>& wfine() const {
    return wfine_;
  }

  void wfine(const blitz::Array<double, 3>& wfine_in) {
    wfine_ = wfine_in;
  }

  

  const blitz::Array<double, 3>& csqfine() const {
    return csqfine_;
  }

  void csqfine(const blitz::Array<double, 3>& csqfine_in) {
    csqfine_ = csqfine_in;
  }

  

  const blitz::Array<double, 3>& cotfine() const {
    return cotfine_;
  }

  void cotfine(const blitz::Array<double, 3>& cotfine_in) {
    cotfine_ = cotfine_in;
  }

  

  const blitz::Array<double, 3>& alphafine() const {
    return alphafine_;
  }

  void alphafine(const blitz::Array<double, 3>& alphafine_in) {
    alphafine_ = alphafine_in;
  }

  

  const blitz::Array<double, 3>& radiifine() const {
    return radiifine_;
  }

  void radiifine(const blitz::Array<double, 3>& radiifine_in) {
    radiifine_ = radiifine_in;
  }

  

  const blitz::Array<int, 1>& ncrit() const {
    return ncrit_;
  }

  

  const blitz::Array<double, 1>& radcrit() const {
    return radcrit_;
  }

  

  const blitz::Array<double, 1>& cotcrit() const {
    return cotcrit_;
  }

  

  const blitz::Array<double, 1>& mu0() const {
    return mu0_;
  }

  

  const blitz::Array<double, 1>& mu1() const {
    return mu1_;
  }

  

  const blitz::Array<double, 1>& cosscat() const {
    return cosscat_;
  }

  

  const blitz::Array<double, 3>& chapfacs() const {
    return chapfacs_;
  }

  

  const blitz::Array<double, 3>& sunpaths() const {
    return sunpaths_;
  }

  

  const blitz::Array<int, 2>& ntraverse() const {
    return ntraverse_;
  }

  

  const blitz::Array<double, 4>& sunpathsfine() const {
    return sunpathsfine_;
  }

  

  const blitz::Array<int, 3>& ntraversefine() const {
    return ntraversefine_;
  }

  

  const bool& fail() const {
    return fail_;
  }

  

  std::string message() const {
    std::string message_ret;
    message_ret = ( std::string(std::string(message_(blitz::Range::all()).begin(), message_(blitz::Range::all()).end()).c_str()) );
    return message_ret;
  }

  

  std::string trace() const {
    std::string trace_ret;
    trace_ret = ( std::string(std::string(trace_(blitz::Range::all()).begin(), trace_(blitz::Range::all()).end()).c_str()) );
    return trace_ret;
  }

  

  
  void run() {
    int message_len = (int) message_.extent(0) - 1;
    int trace_len = (int) trace_.extent(0) - 1;
    
    fo_ssgeometry_master_m_fo_ssgeometry_master_wrap(&maxgeoms_, &maxszas_, &maxvzas_, &maxazms_, &maxlayers_, &maxfine_, &do_obsgeom_, &do_chapman_, &do_planpar_, &do_enhanced_ps_, &ngeoms_, &nszas_, &nvzas_, &nazms_, &nlayers_, &nfine_, &dtr_, &pie_, &vsign_, &eradius_, heights_.dataFirst(), obsgeom_boa_.dataFirst(), alpha_boa_.dataFirst(), theta_boa_.dataFirst(), phi_boa_.dataFirst(), donadir_.dataFirst(), &docrit_, &acrit_, extinc_.dataFirst(), raycon_.dataFirst(), radii_.dataFirst(), alpha_.dataFirst(), cota_.dataFirst(), nfinedivs_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), csqfine_.dataFirst(), cotfine_.dataFirst(), alphafine_.dataFirst(), radiifine_.dataFirst(), ncrit_.dataFirst(), radcrit_.dataFirst(), cotcrit_.dataFirst(), mu0_.dataFirst(), mu1_.dataFirst(), cosscat_.dataFirst(), chapfacs_.dataFirst(), sunpaths_.dataFirst(), ntraverse_.dataFirst(), sunpathsfine_.dataFirst(), ntraversefine_.dataFirst(), &fail_, &message_len, message_.dataFirst(), &trace_len, trace_.dataFirst());
    
  }

  // MANUAL CHANGE
  void print(std::ostream &output_stream) const {
    output_stream << "Fo_Ssgeometry_Master:" << std::endl
      << "      maxgeoms: " << maxgeoms()  << std::endl
      << "       maxszas: " << maxszas()  << std::endl
      << "       maxvzas: " << maxvzas()  << std::endl
      << "       maxazms: " << maxazms()  << std::endl
      << "     maxlayers: " << maxlayers()  << std::endl
      << "       maxfine: " << maxfine()  << std::endl
      << "    do_obsgeom: " << do_obsgeom()  << std::endl
      << "    do_chapman: " << do_chapman()  << std::endl
      << "    do_planpar: " << do_planpar()  << std::endl
      << "do_enhanced_ps: " << do_enhanced_ps()  << std::endl
      << "        ngeoms: " << ngeoms()  << std::endl
      << "         nszas: " << nszas()  << std::endl
      << "         nvzas: " << nvzas()  << std::endl
      << "         nazms: " << nazms()  << std::endl
      << "       nlayers: " << nlayers()  << std::endl
      << "         nfine: " << nfine()  << std::endl
      << "           dtr: " << dtr()  << std::endl
      << "           pie: " << pie()  << std::endl
      << "         vsign: " << vsign()  << std::endl
      << "       eradius: " << eradius()  << std::endl
      << "       heights: " << std::endl << heights()  << std::endl
      << "   obsgeom_boa: " << std::endl << obsgeom_boa()  << std::endl
      << "     alpha_boa: " << std::endl << alpha_boa()  << std::endl
      << "     theta_boa: " << std::endl << theta_boa()  << std::endl
      << "       phi_boa: " << std::endl << phi_boa()  << std::endl
      << "       donadir: " << std::endl << donadir()  << std::endl
      << "        docrit: " << docrit()  << std::endl
      << "         acrit: " << acrit()  << std::endl
      << "        extinc: " << std::endl << extinc()  << std::endl
      << "        raycon: " << std::endl << raycon()  << std::endl
      << "         radii: " << std::endl << radii()  << std::endl
      << "         alpha: " << std::endl << alpha()  << std::endl
      << "          cota: " << std::endl << cota()  << std::endl
      << "     nfinedivs: " << std::endl << nfinedivs()  << std::endl
      << "         xfine: " << std::endl << xfine()  << std::endl
      << "         wfine: " << std::endl << wfine()  << std::endl
      << "       csqfine: " << std::endl << csqfine()  << std::endl
      << "       cotfine: " << std::endl << cotfine()  << std::endl
      << "     alphafine: " << std::endl << alphafine()  << std::endl
      << "     radiifine: " << std::endl << radiifine()  << std::endl
      << "         ncrit: " << std::endl << ncrit()  << std::endl
      << "       radcrit: " << std::endl << radcrit()  << std::endl
      << "       cotcrit: " << std::endl << cotcrit()  << std::endl
      << "           mu0: " << std::endl << mu0()  << std::endl
      << "           mu1: " << std::endl << mu1()  << std::endl
      << "       cosscat: " << std::endl << cosscat()  << std::endl
      << "      chapfacs: " << std::endl << chapfacs()  << std::endl
      << "      sunpaths: " << std::endl << sunpaths()  << std::endl
      << "     ntraverse: " << std::endl << ntraverse()  << std::endl
      << "  sunpathsfine: " << std::endl << sunpathsfine()  << std::endl
      << " ntraversefine: " << std::endl << ntraversefine()  << std::endl
      << "          fail: " << fail()  << std::endl
      << "       message: " << "\"" << message() << "\"" << std::endl
      << "         trace: " << "\"" << trace() << "\"" << std::endl;
  }
  // MANUAL CHANGE

private:
  int maxgeoms_;
  int maxszas_;
  int maxvzas_;
  int maxazms_;
  int maxlayers_;
  int maxfine_;
  bool do_obsgeom_;
  bool do_chapman_;
  bool do_planpar_;
  bool do_enhanced_ps_;
  int ngeoms_;
  int nszas_;
  int nvzas_;
  int nazms_;
  int nlayers_;
  int nfine_;
  double dtr_;
  double pie_;
  double vsign_;
  double eradius_;
  blitz::Array<double, 1> heights_;
  blitz::Array<double, 2> obsgeom_boa_;
  blitz::Array<double, 1> alpha_boa_;
  blitz::Array<double, 1> theta_boa_;
  blitz::Array<double, 1> phi_boa_;
  blitz::Array<bool, 1> donadir_;
  bool docrit_;
  double acrit_;
  blitz::Array<double, 1> extinc_;
  blitz::Array<double, 1> raycon_;
  blitz::Array<double, 1> radii_;
  blitz::Array<double, 2> alpha_;
  blitz::Array<double, 2> cota_;
  blitz::Array<int, 2> nfinedivs_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> csqfine_;
  blitz::Array<double, 3> cotfine_;
  blitz::Array<double, 3> alphafine_;
  blitz::Array<double, 3> radiifine_;
  blitz::Array<int, 1> ncrit_;
  blitz::Array<double, 1> radcrit_;
  blitz::Array<double, 1> cotcrit_;
  blitz::Array<double, 1> mu0_;
  blitz::Array<double, 1> mu1_;
  blitz::Array<double, 1> cosscat_;
  blitz::Array<double, 3> chapfacs_;
  blitz::Array<double, 3> sunpaths_;
  blitz::Array<int, 2> ntraverse_;
  blitz::Array<double, 4> sunpathsfine_;
  blitz::Array<int, 3> ntraversefine_;
  bool fail_;
  blitz::Array<char, 1> message_;
  blitz::Array<char, 1> trace_;
  // MANUAL CHANGE
  Fo_Ssgeometry_Master() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
  // MANUAL CHANGE
};

//-----------------------------------------------------------------------
// Links to module: "fo_scalarss_rtcalcs_i_optimized_m" in file: "FO_ScalarSS_RTCalcs_I.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_scalarss_rtcalcs_i_m_ss_integral_i_up_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxfine_in, const bool* do_planpar_in, const bool* do_regular_ps_in, const bool* do_enhanced_ps_in, const bool* donadir_in, const bool* do_sleave_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* aclevel_in, const double* reflec_in, const double* slterm_in, const double* extinction_in, const double* deltaus_in, const double* exactscat_up_in, const double* flux_in, const double* mu0_in, const double* mu1_in, const int* ncrit_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* raycon_in, const double* cota_in, const double* sunpaths_in, const int* ntraverse_in, const double* sunpathsfine_in, const int* ntraversefine_in, const double* intensity_up_in, const double* intensity_db_in, const double* cumsource_up_in);
}

class Fo_Scalarss_Rtcalcs_I : public virtual GenericObject {

public:
  Fo_Scalarss_Rtcalcs_I(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfine_in, const int& ngeoms_in, const int& nlayers_in) : maxgeoms_(maxgeoms_in), maxlayers_(maxlayers_in), maxfine_(maxfine_in), ngeoms_(ngeoms_in), nlayers_(nlayers_in) 
  { 
    do_planpar_ = false;
    do_regular_ps_ = false;
    do_enhanced_ps_ = false;
    donadir_.reference( blitz::Array<bool, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    donadir_ = false;
    do_sleave_ = false;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    aclevel_ = 0;
    reflec_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    reflec_ = 0;
    slterm_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    slterm_ = 0;
    extinction_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    extinction_ = 0;
    deltaus_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    deltaus_ = 0;
    exactscat_up_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    exactscat_up_ = 0;
    flux_ = 0;
    mu0_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu0_ = 0;
    mu1_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu1_ = 0;
    ncrit_.reference( blitz::Array<int, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    ncrit_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    csqfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    csqfine_ = 0;
    cotfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    cotfine_ = 0;
    raycon_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    raycon_ = 0;
    cota_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cota_ = 0;
    sunpaths_.reference( blitz::Array<double, 3>(maxlayers_-0+1, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_ = 0;
    ntraverse_.reference( blitz::Array<int, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_ = 0;
    sunpathsfine_.reference( blitz::Array<double, 4>(maxlayers_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_ = 0;
    ntraversefine_.reference( blitz::Array<int, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_ = 0;
    intensity_up_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    intensity_up_ = 0;
    intensity_db_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    intensity_db_ = 0;
    cumsource_up_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cumsource_up_ = 0;
    // Initialize type pointers
    
  }

  virtual ~Fo_Scalarss_Rtcalcs_I() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxfine() const {
    return maxfine_;
  }

  

  const bool& do_planpar() const {
    return do_planpar_;
  }

  void do_planpar(const bool& do_planpar_in) {
    do_planpar_ = do_planpar_in;
  }

  

  const bool& do_regular_ps() const {
    return do_regular_ps_;
  }

  void do_regular_ps(const bool& do_regular_ps_in) {
    do_regular_ps_ = do_regular_ps_in;
  }

  

  const bool& do_enhanced_ps() const {
    return do_enhanced_ps_;
  }

  void do_enhanced_ps(const bool& do_enhanced_ps_in) {
    do_enhanced_ps_ = do_enhanced_ps_in;
  }

  

  const blitz::Array<bool, 1>& donadir() const {
    return donadir_;
  }

  void donadir(const blitz::Array<bool, 1>& donadir_in) {
    donadir_ = donadir_in;
  }

  

  const bool& do_sleave() const {
    return do_sleave_;
  }

  void do_sleave(const bool& do_sleave_in) {
    do_sleave_ = do_sleave_in;
  }

  

  const int& ngeoms() const {
    return ngeoms_;
  }

  

  const int& nlayers() const {
    return nlayers_;
  }

  

  const blitz::Array<int, 2>& nfinedivs() const {
    return nfinedivs_;
  }

  void nfinedivs(const blitz::Array<int, 2>& nfinedivs_in) {
    nfinedivs_ = nfinedivs_in;
  }

  

  const int& aclevel() const {
    return aclevel_;
  }

  void aclevel(const int& aclevel_in) {
    aclevel_ = aclevel_in;
  }

  

  const blitz::Array<double, 1>& reflec() const {
    return reflec_;
  }

  void reflec(const blitz::Array<double, 1>& reflec_in) {
    reflec_ = reflec_in;
  }

  

  const blitz::Array<double, 1>& slterm() const {
    return slterm_;
  }

  void slterm(const blitz::Array<double, 1>& slterm_in) {
    slterm_ = slterm_in;
  }

  

  const blitz::Array<double, 1>& extinction() const {
    return extinction_;
  }

  void extinction(const blitz::Array<double, 1>& extinction_in) {
    extinction_ = extinction_in;
  }

  

  const blitz::Array<double, 1>& deltaus() const {
    return deltaus_;
  }

  void deltaus(const blitz::Array<double, 1>& deltaus_in) {
    deltaus_ = deltaus_in;
  }

  

  const blitz::Array<double, 2>& exactscat_up() const {
    return exactscat_up_;
  }

  void exactscat_up(const blitz::Array<double, 2>& exactscat_up_in) {
    exactscat_up_ = exactscat_up_in;
  }

  

  const double& flux() const {
    return flux_;
  }

  void flux(const double& flux_in) {
    flux_ = flux_in;
  }

  

  const blitz::Array<double, 1>& mu0() const {
    return mu0_;
  }

  void mu0(const blitz::Array<double, 1>& mu0_in) {
    mu0_ = mu0_in;
  }

  

  const blitz::Array<double, 1>& mu1() const {
    return mu1_;
  }

  void mu1(const blitz::Array<double, 1>& mu1_in) {
    mu1_ = mu1_in;
  }

  

  const blitz::Array<int, 1>& ncrit() const {
    return ncrit_;
  }

  void ncrit(const blitz::Array<int, 1>& ncrit_in) {
    ncrit_ = ncrit_in;
  }

  

  const blitz::Array<double, 3>& xfine() const {
    return xfine_;
  }

  void xfine(const blitz::Array<double, 3>& xfine_in) {
    xfine_ = xfine_in;
  }

  

  const blitz::Array<double, 3>& wfine() const {
    return wfine_;
  }

  void wfine(const blitz::Array<double, 3>& wfine_in) {
    wfine_ = wfine_in;
  }

  

  const blitz::Array<double, 3>& csqfine() const {
    return csqfine_;
  }

  void csqfine(const blitz::Array<double, 3>& csqfine_in) {
    csqfine_ = csqfine_in;
  }

  

  const blitz::Array<double, 3>& cotfine() const {
    return cotfine_;
  }

  void cotfine(const blitz::Array<double, 3>& cotfine_in) {
    cotfine_ = cotfine_in;
  }

  

  const blitz::Array<double, 1>& raycon() const {
    return raycon_;
  }

  void raycon(const blitz::Array<double, 1>& raycon_in) {
    raycon_ = raycon_in;
  }

  

  const blitz::Array<double, 2>& cota() const {
    return cota_;
  }

  void cota(const blitz::Array<double, 2>& cota_in) {
    cota_ = cota_in;
  }

  

  const blitz::Array<double, 3>& sunpaths() const {
    return sunpaths_;
  }

  void sunpaths(const blitz::Array<double, 3>& sunpaths_in) {
    sunpaths_ = sunpaths_in;
  }

  

  const blitz::Array<int, 2>& ntraverse() const {
    return ntraverse_;
  }

  void ntraverse(const blitz::Array<int, 2>& ntraverse_in) {
    ntraverse_ = ntraverse_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine() const {
    return sunpathsfine_;
  }

  void sunpathsfine(const blitz::Array<double, 4>& sunpathsfine_in) {
    sunpathsfine_ = sunpathsfine_in;
  }

  

  const blitz::Array<int, 3>& ntraversefine() const {
    return ntraversefine_;
  }

  void ntraversefine(const blitz::Array<int, 3>& ntraversefine_in) {
    ntraversefine_ = ntraversefine_in;
  }

  

  const blitz::Array<double, 1>& intensity_up() const {
    return intensity_up_;
  }

  

  const blitz::Array<double, 1>& intensity_db() const {
    return intensity_db_;
  }

  

  const blitz::Array<double, 2>& cumsource_up() const {
    return cumsource_up_;
  }

  

  
  void ss_integral_i_up() {
    
    
    fo_scalarss_rtcalcs_i_m_ss_integral_i_up_wrap(&maxgeoms_, &maxlayers_, &maxfine_, &do_planpar_, &do_regular_ps_, &do_enhanced_ps_, donadir_.dataFirst(), &do_sleave_, &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &aclevel_, reflec_.dataFirst(), slterm_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), exactscat_up_.dataFirst(), &flux_, mu0_.dataFirst(), mu1_.dataFirst(), ncrit_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), csqfine_.dataFirst(), cotfine_.dataFirst(), raycon_.dataFirst(), cota_.dataFirst(), sunpaths_.dataFirst(), ntraverse_.dataFirst(), sunpathsfine_.dataFirst(), ntraversefine_.dataFirst(), intensity_up_.dataFirst(), intensity_db_.dataFirst(), cumsource_up_.dataFirst());
    
  }

  friend std::ostream& operator<<(std::ostream &output_stream, const Fo_Scalarss_Rtcalcs_I &obj) {
    output_stream << "Fo_Scalarss_Rtcalcs_I:" << std::endl
      << "      maxgeoms: " << obj.maxgeoms()  << std::endl
      << "     maxlayers: " << obj.maxlayers()  << std::endl
      << "       maxfine: " << obj.maxfine()  << std::endl
      << "    do_planpar: " << obj.do_planpar()  << std::endl
      << " do_regular_ps: " << obj.do_regular_ps()  << std::endl
      << "do_enhanced_ps: " << obj.do_enhanced_ps()  << std::endl
      << "       donadir: " << std::endl << obj.donadir()  << std::endl
      << "     do_sleave: " << obj.do_sleave()  << std::endl
      << "        ngeoms: " << obj.ngeoms()  << std::endl
      << "       nlayers: " << obj.nlayers()  << std::endl
      << "     nfinedivs: " << std::endl << obj.nfinedivs()  << std::endl
      << "       aclevel: " << obj.aclevel()  << std::endl
      << "        reflec: " << std::endl << obj.reflec()  << std::endl
      << "        slterm: " << std::endl << obj.slterm()  << std::endl
      << "    extinction: " << std::endl << obj.extinction()  << std::endl
      << "       deltaus: " << std::endl << obj.deltaus()  << std::endl
      << "  exactscat_up: " << std::endl << obj.exactscat_up()  << std::endl
      << "          flux: " << obj.flux()  << std::endl
      << "           mu0: " << std::endl << obj.mu0()  << std::endl
      << "           mu1: " << std::endl << obj.mu1()  << std::endl
      << "         ncrit: " << std::endl << obj.ncrit()  << std::endl
      << "         xfine: " << std::endl << obj.xfine()  << std::endl
      << "         wfine: " << std::endl << obj.wfine()  << std::endl
      << "       csqfine: " << std::endl << obj.csqfine()  << std::endl
      << "       cotfine: " << std::endl << obj.cotfine()  << std::endl
      << "        raycon: " << std::endl << obj.raycon()  << std::endl
      << "          cota: " << std::endl << obj.cota()  << std::endl
      << "      sunpaths: " << std::endl << obj.sunpaths()  << std::endl
      << "     ntraverse: " << std::endl << obj.ntraverse()  << std::endl
      << "  sunpathsfine: " << std::endl << obj.sunpathsfine()  << std::endl
      << " ntraversefine: " << std::endl << obj.ntraversefine()  << std::endl
      << "  intensity_up: " << std::endl << obj.intensity_up()  << std::endl
      << "  intensity_db: " << std::endl << obj.intensity_db()  << std::endl
      << "  cumsource_up: " << std::endl << obj.cumsource_up()  << std::endl;
    return output_stream;

  }

private:
  int maxgeoms_;
  int maxlayers_;
  int maxfine_;
  bool do_planpar_;
  bool do_regular_ps_;
  bool do_enhanced_ps_;
  blitz::Array<bool, 1> donadir_;
  bool do_sleave_;
  int ngeoms_;
  int nlayers_;
  blitz::Array<int, 2> nfinedivs_;
  int aclevel_;
  blitz::Array<double, 1> reflec_;
  blitz::Array<double, 1> slterm_;
  blitz::Array<double, 1> extinction_;
  blitz::Array<double, 1> deltaus_;
  blitz::Array<double, 2> exactscat_up_;
  double flux_;
  blitz::Array<double, 1> mu0_;
  blitz::Array<double, 1> mu1_;
  blitz::Array<int, 1> ncrit_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> csqfine_;
  blitz::Array<double, 3> cotfine_;
  blitz::Array<double, 1> raycon_;
  blitz::Array<double, 2> cota_;
  blitz::Array<double, 3> sunpaths_;
  blitz::Array<int, 2> ntraverse_;
  blitz::Array<double, 4> sunpathsfine_;
  blitz::Array<int, 3> ntraversefine_;
  blitz::Array<double, 1> intensity_up_;
  blitz::Array<double, 1> intensity_db_;
  blitz::Array<double, 2> cumsource_up_;
};

//-----------------------------------------------------------------------
// Links to module: "fo_scalarss_rtcalcs_ilps_optimized_m" in file: "FO_ScalarSS_RTCalcs_ILPS.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_up_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxfine_in, const int* max_atmoswfs_in, const int* max_surfacewfs_in, const bool* do_planpar_in, const bool* do_regular_ps_in, const bool* do_enhanced_ps_in, const bool* donadir_in, const bool* do_sleave_in, const bool* do_profilewfs_in, const bool* do_reflecwfs_in, const bool* do_sleavewfs_in, const bool* lvaryflags_in, const int* lvarynums_in, const int* n_reflecwfs_in, const int* n_sleavewfs_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* aclevel_in, const double* reflec_in, const double* slterm_in, const double* extinction_in, const double* deltaus_in, const double* exactscat_up_in, const double* flux_in, const double* ls_reflec_in, const double* lssl_slterm_in, const double* l_extinction_in, const double* l_deltaus_in, const double* l_exactscat_up_in, const double* mu0_in, const double* mu1_in, const int* ncrit_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* raycon_in, const double* cota_in, const double* sunpaths_in, const int* ntraverse_in, const double* sunpaths_fine_in, const int* ntraverse_fine_in, const double* intensity_up_in, const double* intensity_db_in, const double* lp_jacobians_up_in, const double* lp_jacobians_db_in, const double* ls_jacobians_db_in);

  void first_order_debug_output(const int* maxgeoms_in, const int* maxlayers_in, const int* maxfine_in, const bool* do_planpar_in, const bool* do_regular_ps_in, const bool* do_enhanced_ps_in, const bool* donadir_in, const bool* do_sleave_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* aclevel_in, const double* reflec_in, const double* slterm_in, const double* extinction_in, const double* deltaus_in, const double* exactscat_up_in, const double* flux_in, const double* mu0_in, const double* mu1_in, const int* ncrit_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* raycon_in, const double* cota_in, const double* sunpaths_in, const int* ntraverse_in, const double* sunpathsfine_in, const int* ntraversefine_in);
}

  // MANUAL CHANGE
class Fo_Scalarss_Rtcalcs_Ilps : public Printable<Fo_Scalarss_Rtcalcs_Ilps> {
  // MANUAL CHANGE

public:
  Fo_Scalarss_Rtcalcs_Ilps(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfine_in, const int& max_atmoswfs_in, const int& max_surfacewfs_in, const int& ngeoms_in, const int& nlayers_in) : maxgeoms_(maxgeoms_in), maxlayers_(maxlayers_in), maxfine_(maxfine_in), max_atmoswfs_(max_atmoswfs_in), max_surfacewfs_(max_surfacewfs_in), ngeoms_(ngeoms_in), nlayers_(nlayers_in) 
  { 
    do_planpar_ = false;
    do_regular_ps_ = false;
    do_enhanced_ps_ = false;
    donadir_.reference( blitz::Array<bool, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    donadir_ = false;
    do_sleave_ = false;
    do_profilewfs_ = false;
    do_reflecwfs_ = false;
    do_sleavewfs_ = false;
    lvaryflags_.reference( blitz::Array<bool, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    lvaryflags_ = false;
    lvarynums_.reference( blitz::Array<int, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    lvarynums_ = 0;
    n_reflecwfs_ = 0;
    n_sleavewfs_ = 0;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    aclevel_ = 0;
    reflec_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    reflec_ = 0;
    slterm_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    slterm_ = 0;
    extinction_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    extinction_ = 0;
    deltaus_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    deltaus_ = 0;
    exactscat_up_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    exactscat_up_ = 0;
    flux_ = 0;
    ls_reflec_.reference( blitz::Array<double, 2>(maxgeoms_, max_surfacewfs_, blitz::ColumnMajorArray<2>()) );
    ls_reflec_ = 0;
    lssl_slterm_.reference( blitz::Array<double, 2>(maxgeoms_, max_surfacewfs_, blitz::ColumnMajorArray<2>()) );
    lssl_slterm_ = 0;
    l_extinction_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_extinction_ = 0;
    l_deltaus_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_deltaus_ = 0;
    l_exactscat_up_.reference( blitz::Array<double, 3>(maxlayers_, maxgeoms_, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    l_exactscat_up_ = 0;
    mu0_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu0_ = 0;
    mu1_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu1_ = 0;
    ncrit_.reference( blitz::Array<int, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    ncrit_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    csqfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    csqfine_ = 0;
    cotfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    cotfine_ = 0;
    raycon_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    raycon_ = 0;
    cota_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cota_ = 0;
    sunpaths_.reference( blitz::Array<double, 3>(maxlayers_-0+1, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_ = 0;
    ntraverse_.reference( blitz::Array<int, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_ = 0;
    sunpaths_fine_.reference( blitz::Array<double, 4>(maxlayers_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpaths_fine_ = 0;
    ntraverse_fine_.reference( blitz::Array<int, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraverse_fine_ = 0;
    intensity_up_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    intensity_up_ = 0;
    intensity_db_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    intensity_db_ = 0;
    lp_jacobians_up_.reference( blitz::Array<double, 3>(maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    lp_jacobians_up_ = 0;
    lp_jacobians_db_.reference( blitz::Array<double, 3>(maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    lp_jacobians_db_ = 0;
    ls_jacobians_db_.reference( blitz::Array<double, 2>(maxgeoms_, max_surfacewfs_, blitz::ColumnMajorArray<2>()) );
    ls_jacobians_db_ = 0;
    // Initialize type pointers
    
  }

  virtual ~Fo_Scalarss_Rtcalcs_Ilps() = default;

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxfine() const {
    return maxfine_;
  }

  

  const int& max_atmoswfs() const {
    return max_atmoswfs_;
  }

  

  const int& max_surfacewfs() const {
    return max_surfacewfs_;
  }

  

  const bool& do_planpar() const {
    return do_planpar_;
  }

  void do_planpar(const bool& do_planpar_in) {
    do_planpar_ = do_planpar_in;
  }

  

  const bool& do_regular_ps() const {
    return do_regular_ps_;
  }

  void do_regular_ps(const bool& do_regular_ps_in) {
    do_regular_ps_ = do_regular_ps_in;
  }

  

  const bool& do_enhanced_ps() const {
    return do_enhanced_ps_;
  }

  void do_enhanced_ps(const bool& do_enhanced_ps_in) {
    do_enhanced_ps_ = do_enhanced_ps_in;
  }

  

  const blitz::Array<bool, 1>& donadir() const {
    return donadir_;
  }

  void donadir(const blitz::Array<bool, 1>& donadir_in) {
    donadir_ = donadir_in;
  }

  

  const bool& do_sleave() const {
    return do_sleave_;
  }

  void do_sleave(const bool& do_sleave_in) {
    do_sleave_ = do_sleave_in;
  }

  

  const bool& do_profilewfs() const {
    return do_profilewfs_;
  }

  void do_profilewfs(const bool& do_profilewfs_in) {
    do_profilewfs_ = do_profilewfs_in;
  }

  

  const bool& do_reflecwfs() const {
    return do_reflecwfs_;
  }

  void do_reflecwfs(const bool& do_reflecwfs_in) {
    do_reflecwfs_ = do_reflecwfs_in;
  }

  

  const bool& do_sleavewfs() const {
    return do_sleavewfs_;
  }

  void do_sleavewfs(const bool& do_sleavewfs_in) {
    do_sleavewfs_ = do_sleavewfs_in;
  }

  

  const blitz::Array<bool, 1>& lvaryflags() const {
    return lvaryflags_;
  }

  void lvaryflags(const blitz::Array<bool, 1>& lvaryflags_in) {
    lvaryflags_ = lvaryflags_in;
  }

  

  const blitz::Array<int, 1>& lvarynums() const {
    return lvarynums_;
  }

  void lvarynums(const blitz::Array<int, 1>& lvarynums_in) {
    lvarynums_ = lvarynums_in;
  }

  

  const int& n_reflecwfs() const {
    return n_reflecwfs_;
  }

  void n_reflecwfs(const int& n_reflecwfs_in) {
    n_reflecwfs_ = n_reflecwfs_in;
  }

  

  const int& n_sleavewfs() const {
    return n_sleavewfs_;
  }

  void n_sleavewfs(const int& n_sleavewfs_in) {
    n_sleavewfs_ = n_sleavewfs_in;
  }

  

  const int& ngeoms() const {
    return ngeoms_;
  }

  

  const int& nlayers() const {
    return nlayers_;
  }

  

  const blitz::Array<int, 2>& nfinedivs() const {
    return nfinedivs_;
  }

  void nfinedivs(const blitz::Array<int, 2>& nfinedivs_in) {
    nfinedivs_ = nfinedivs_in;
  }

  

  const int& aclevel() const {
    return aclevel_;
  }

  void aclevel(const int& aclevel_in) {
    aclevel_ = aclevel_in;
  }

  

  const blitz::Array<double, 1>& reflec() const {
    return reflec_;
  }

  void reflec(const blitz::Array<double, 1>& reflec_in) {
    reflec_ = reflec_in;
  }

  

  const blitz::Array<double, 1>& slterm() const {
    return slterm_;
  }

  void slterm(const blitz::Array<double, 1>& slterm_in) {
    slterm_ = slterm_in;
  }

  

  const blitz::Array<double, 1>& extinction() const {
    return extinction_;
  }

  void extinction(const blitz::Array<double, 1>& extinction_in) {
    extinction_ = extinction_in;
  }

  

  const blitz::Array<double, 1>& deltaus() const {
    return deltaus_;
  }

  void deltaus(const blitz::Array<double, 1>& deltaus_in) {
    deltaus_ = deltaus_in;
  }

  

  const blitz::Array<double, 2>& exactscat_up() const {
    return exactscat_up_;
  }

  void exactscat_up(const blitz::Array<double, 2>& exactscat_up_in) {
    exactscat_up_ = exactscat_up_in;
  }

  

  const double& flux() const {
    return flux_;
  }

  void flux(const double& flux_in) {
    flux_ = flux_in;
  }

  

  const blitz::Array<double, 2>& ls_reflec() const {
    return ls_reflec_;
  }

  void ls_reflec(const blitz::Array<double, 2>& ls_reflec_in) {
    ls_reflec_ = ls_reflec_in;
  }

  

  const blitz::Array<double, 2>& lssl_slterm() const {
    return lssl_slterm_;
  }

  void lssl_slterm(const blitz::Array<double, 2>& lssl_slterm_in) {
    lssl_slterm_ = lssl_slterm_in;
  }

  

  const blitz::Array<double, 2>& l_extinction() const {
    return l_extinction_;
  }

  void l_extinction(const blitz::Array<double, 2>& l_extinction_in) {
    l_extinction_ = l_extinction_in;
  }

  

  const blitz::Array<double, 2>& l_deltaus() const {
    return l_deltaus_;
  }

  void l_deltaus(const blitz::Array<double, 2>& l_deltaus_in) {
    l_deltaus_ = l_deltaus_in;
  }

  

  const blitz::Array<double, 3>& l_exactscat_up() const {
    return l_exactscat_up_;
  }

  void l_exactscat_up(const blitz::Array<double, 3>& l_exactscat_up_in) {
    l_exactscat_up_ = l_exactscat_up_in;
  }

  

  const blitz::Array<double, 1>& mu0() const {
    return mu0_;
  }

  void mu0(const blitz::Array<double, 1>& mu0_in) {
    mu0_ = mu0_in;
  }

  

  const blitz::Array<double, 1>& mu1() const {
    return mu1_;
  }

  void mu1(const blitz::Array<double, 1>& mu1_in) {
    mu1_ = mu1_in;
  }

  

  const blitz::Array<int, 1>& ncrit() const {
    return ncrit_;
  }

  void ncrit(const blitz::Array<int, 1>& ncrit_in) {
    ncrit_ = ncrit_in;
  }

  

  const blitz::Array<double, 3>& xfine() const {
    return xfine_;
  }

  void xfine(const blitz::Array<double, 3>& xfine_in) {
    xfine_ = xfine_in;
  }

  

  const blitz::Array<double, 3>& wfine() const {
    return wfine_;
  }

  void wfine(const blitz::Array<double, 3>& wfine_in) {
    wfine_ = wfine_in;
  }

  

  const blitz::Array<double, 3>& csqfine() const {
    return csqfine_;
  }

  void csqfine(const blitz::Array<double, 3>& csqfine_in) {
    csqfine_ = csqfine_in;
  }

  

  const blitz::Array<double, 3>& cotfine() const {
    return cotfine_;
  }

  void cotfine(const blitz::Array<double, 3>& cotfine_in) {
    cotfine_ = cotfine_in;
  }

  

  const blitz::Array<double, 1>& raycon() const {
    return raycon_;
  }

  void raycon(const blitz::Array<double, 1>& raycon_in) {
    raycon_ = raycon_in;
  }

  

  const blitz::Array<double, 2>& cota() const {
    return cota_;
  }

  void cota(const blitz::Array<double, 2>& cota_in) {
    cota_ = cota_in;
  }

  

  const blitz::Array<double, 3>& sunpaths() const {
    return sunpaths_;
  }

  void sunpaths(const blitz::Array<double, 3>& sunpaths_in) {
    sunpaths_ = sunpaths_in;
  }

  

  const blitz::Array<int, 2>& ntraverse() const {
    return ntraverse_;
  }

  void ntraverse(const blitz::Array<int, 2>& ntraverse_in) {
    ntraverse_ = ntraverse_in;
  }

  

  const blitz::Array<double, 4>& sunpaths_fine() const {
    return sunpaths_fine_;
  }

  void sunpaths_fine(const blitz::Array<double, 4>& sunpaths_fine_in) {
    sunpaths_fine_ = sunpaths_fine_in;
  }

  

  const blitz::Array<int, 3>& ntraverse_fine() const {
    return ntraverse_fine_;
  }

  void ntraverse_fine(const blitz::Array<int, 3>& ntraverse_fine_in) {
    ntraverse_fine_ = ntraverse_fine_in;
  }

  

  const blitz::Array<double, 1>& intensity_up() const {
    return intensity_up_;
  }

  

  const blitz::Array<double, 1>& intensity_db() const {
    return intensity_db_;
  }

  

  const blitz::Array<double, 3>& lp_jacobians_up() const {
    return lp_jacobians_up_;
  }

  

  const blitz::Array<double, 3>& lp_jacobians_db() const {
    return lp_jacobians_db_;
  }

  

  const blitz::Array<double, 2>& ls_jacobians_db() const {
    return ls_jacobians_db_;
  }

  

  
  void ss_integral_ilps_up() {

#if FIRST_ORDER_DEBUG
    first_order_debug_output(&maxgeoms_, &maxlayers_, &maxfine_, &do_planpar_, &do_regular_ps_, &do_enhanced_ps_, donadir_.dataFirst(), &do_sleave_, &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &aclevel_, reflec_.dataFirst(), slterm_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), exactscat_up_.dataFirst(), &flux_, mu0_.dataFirst(), mu1_.dataFirst(), ncrit_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), csqfine_.dataFirst(), cotfine_.dataFirst(), raycon_.dataFirst(), cota_.dataFirst(), sunpaths_.dataFirst(), ntraverse_.dataFirst(), sunpaths_fine_.dataFirst(), ntraverse_fine_.dataFirst());
#endif
    
    fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_up_wrap(&maxgeoms_, &maxlayers_, &maxfine_, &max_atmoswfs_, &max_surfacewfs_, &do_planpar_, &do_regular_ps_, &do_enhanced_ps_, donadir_.dataFirst(), &do_sleave_, &do_profilewfs_, &do_reflecwfs_, &do_sleavewfs_, lvaryflags_.dataFirst(), lvarynums_.dataFirst(), &n_reflecwfs_, &n_sleavewfs_, &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &aclevel_, reflec_.dataFirst(), slterm_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), exactscat_up_.dataFirst(), &flux_, ls_reflec_.dataFirst(), lssl_slterm_.dataFirst(), l_extinction_.dataFirst(), l_deltaus_.dataFirst(), l_exactscat_up_.dataFirst(), mu0_.dataFirst(), mu1_.dataFirst(), ncrit_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), csqfine_.dataFirst(), cotfine_.dataFirst(), raycon_.dataFirst(), cota_.dataFirst(), sunpaths_.dataFirst(), ntraverse_.dataFirst(), sunpaths_fine_.dataFirst(), ntraverse_fine_.dataFirst(), intensity_up_.dataFirst(), intensity_db_.dataFirst(), lp_jacobians_up_.dataFirst(), lp_jacobians_db_.dataFirst(), ls_jacobians_db_.dataFirst());
    
  }

  // MANUAL CHANGE
  void print(std::ostream &output_stream) const {
    output_stream << "Fo_Scalarss_Rtcalcs_Ilps:" << std::endl
      << "       maxgeoms: " << maxgeoms()  << std::endl
      << "      maxlayers: " << maxlayers()  << std::endl
      << "        maxfine: " << maxfine()  << std::endl
      << "   max_atmoswfs: " << max_atmoswfs()  << std::endl
      << " max_surfacewfs: " << max_surfacewfs()  << std::endl
      << "     do_planpar: " << do_planpar()  << std::endl
      << "  do_regular_ps: " << do_regular_ps()  << std::endl
      << " do_enhanced_ps: " << do_enhanced_ps()  << std::endl
      << "        donadir: " << std::endl << donadir()  << std::endl
      << "      do_sleave: " << do_sleave()  << std::endl
      << "  do_profilewfs: " << do_profilewfs()  << std::endl
      << "   do_reflecwfs: " << do_reflecwfs()  << std::endl
      << "   do_sleavewfs: " << do_sleavewfs()  << std::endl
      << "     lvaryflags: " << std::endl << lvaryflags()  << std::endl
      << "      lvarynums: " << std::endl << lvarynums()  << std::endl
      << "    n_reflecwfs: " << n_reflecwfs()  << std::endl
      << "    n_sleavewfs: " << n_sleavewfs()  << std::endl
      << "         ngeoms: " << ngeoms()  << std::endl
      << "        nlayers: " << nlayers()  << std::endl
      << "      nfinedivs: " << std::endl << nfinedivs()  << std::endl
      << "        aclevel: " << aclevel()  << std::endl
      << "         reflec: " << std::endl << reflec()  << std::endl
      << "         slterm: " << std::endl << slterm()  << std::endl
      << "     extinction: " << std::endl << extinction()  << std::endl
      << "        deltaus: " << std::endl << deltaus()  << std::endl
      << "   exactscat_up: " << std::endl << exactscat_up()  << std::endl
      << "           flux: " << flux()  << std::endl
      << "      ls_reflec: " << std::endl << ls_reflec()  << std::endl
      << "    lssl_slterm: " << std::endl << lssl_slterm()  << std::endl
      << "   l_extinction: " << std::endl << l_extinction()  << std::endl
      << "      l_deltaus: " << std::endl << l_deltaus()  << std::endl
      << " l_exactscat_up: " << std::endl << l_exactscat_up()  << std::endl
      << "            mu0: " << std::endl << mu0()  << std::endl
      << "            mu1: " << std::endl << mu1()  << std::endl
      << "          ncrit: " << std::endl << ncrit()  << std::endl
      << "          xfine: " << std::endl << xfine()  << std::endl
      << "          wfine: " << std::endl << wfine()  << std::endl
      << "        csqfine: " << std::endl << csqfine()  << std::endl
      << "        cotfine: " << std::endl << cotfine()  << std::endl
      << "         raycon: " << std::endl << raycon()  << std::endl
      << "           cota: " << std::endl << cota()  << std::endl
      << "       sunpaths: " << std::endl << sunpaths()  << std::endl
      << "      ntraverse: " << std::endl << ntraverse()  << std::endl
      << "  sunpaths_fine: " << std::endl << sunpaths_fine()  << std::endl
      << " ntraverse_fine: " << std::endl << ntraverse_fine()  << std::endl
      << "   intensity_up: " << std::endl << intensity_up()  << std::endl
      << "   intensity_db: " << std::endl << intensity_db()  << std::endl
      << "lp_jacobians_up: " << std::endl << lp_jacobians_up()  << std::endl
      << "lp_jacobians_db: " << std::endl << lp_jacobians_db()  << std::endl
      << "ls_jacobians_db: " << std::endl << ls_jacobians_db()  << std::endl;
  }
  // MANUAL CHANGE

private:
  int maxgeoms_;
  int maxlayers_;
  int maxfine_;
  int max_atmoswfs_;
  int max_surfacewfs_;
  bool do_planpar_;
  bool do_regular_ps_;
  bool do_enhanced_ps_;
  blitz::Array<bool, 1> donadir_;
  bool do_sleave_;
  bool do_profilewfs_;
  bool do_reflecwfs_;
  bool do_sleavewfs_;
  blitz::Array<bool, 1> lvaryflags_;
  blitz::Array<int, 1> lvarynums_;
  int n_reflecwfs_;
  int n_sleavewfs_;
  int ngeoms_;
  int nlayers_;
  blitz::Array<int, 2> nfinedivs_;
  int aclevel_;
  blitz::Array<double, 1> reflec_;
  blitz::Array<double, 1> slterm_;
  blitz::Array<double, 1> extinction_;
  blitz::Array<double, 1> deltaus_;
  blitz::Array<double, 2> exactscat_up_;
  double flux_;
  blitz::Array<double, 2> ls_reflec_;
  blitz::Array<double, 2> lssl_slterm_;
  blitz::Array<double, 2> l_extinction_;
  blitz::Array<double, 2> l_deltaus_;
  blitz::Array<double, 3> l_exactscat_up_;
  blitz::Array<double, 1> mu0_;
  blitz::Array<double, 1> mu1_;
  blitz::Array<int, 1> ncrit_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> csqfine_;
  blitz::Array<double, 3> cotfine_;
  blitz::Array<double, 1> raycon_;
  blitz::Array<double, 2> cota_;
  blitz::Array<double, 3> sunpaths_;
  blitz::Array<int, 2> ntraverse_;
  blitz::Array<double, 4> sunpaths_fine_;
  blitz::Array<int, 3> ntraverse_fine_;
  blitz::Array<double, 1> intensity_up_;
  blitz::Array<double, 1> intensity_db_;
  blitz::Array<double, 3> lp_jacobians_up_;
  blitz::Array<double, 3> lp_jacobians_db_;
  blitz::Array<double, 2> ls_jacobians_db_;
  // MANUAL CHANGE
  Fo_Scalarss_Rtcalcs_Ilps() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
  // MANUAL CHANGE
};

//-----------------------------------------------------------------------
// Links to module: "fo_scalarss_spherfuncs_m" in file: "FO_ScalarSS_Spherfuncs.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_scalarss_spherfuncs_m_fo_scalarss_spherfuncs_wrap(const bool* starter_in, const int* maxmoms_in, const int* maxgeoms_in, const int* nmoms_in, const int* ngeoms_in, const double* df1_in, const double* df2_in, const double* cosscat_in, const double* ss_pleg_in);
}

  // MANUAL CHANGE
class Fo_Scalarss_Spherfuncs : public Printable<Fo_Scalarss_Spherfuncs> {
  // MANUAL CHANGE

public:
  Fo_Scalarss_Spherfuncs(const bool& starter_in, const int& maxmoms_in, const int& maxgeoms_in, const int& nmoms_in, const int& ngeoms_in) : starter_(starter_in), maxmoms_(maxmoms_in), maxgeoms_(maxgeoms_in), nmoms_(nmoms_in), ngeoms_(ngeoms_in) 
  { 
    df1_.reference( blitz::Array<double, 1>(maxmoms_, blitz::ColumnMajorArray<1>()) );
    df1_ = 0;
    df2_.reference( blitz::Array<double, 1>(maxmoms_, blitz::ColumnMajorArray<1>()) );
    df2_ = 0;
    cosscat_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    cosscat_ = 0;
    ss_pleg_.reference( blitz::Array<double, 2>(maxmoms_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ss_pleg_ = 0;
    // Initialize type pointers
    
  }

  virtual ~Fo_Scalarss_Spherfuncs() = default;

  const bool& starter() const {
    return starter_;
  }

  

  const int& maxmoms() const {
    return maxmoms_;
  }

  

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& nmoms() const {
    return nmoms_;
  }

  

  const int& ngeoms() const {
    return ngeoms_;
  }

  

  const blitz::Array<double, 1>& df1() const {
    return df1_;
  }

  void df1(const blitz::Array<double, 1>& df1_in) {
    df1_ = df1_in;
  }

  

  const blitz::Array<double, 1>& df2() const {
    return df2_;
  }

  void df2(const blitz::Array<double, 1>& df2_in) {
    df2_ = df2_in;
  }

  

  const blitz::Array<double, 1>& cosscat() const {
    return cosscat_;
  }

  void cosscat(const blitz::Array<double, 1>& cosscat_in) {
    cosscat_ = cosscat_in;
  }

  

  const blitz::Array<double, 2>& ss_pleg() const {
    return ss_pleg_;
  }

  

  
  void run() {
    
    
    fo_scalarss_spherfuncs_m_fo_scalarss_spherfuncs_wrap(&starter_, &maxmoms_, &maxgeoms_, &nmoms_, &ngeoms_, df1_.dataFirst(), df2_.dataFirst(), cosscat_.dataFirst(), ss_pleg_.dataFirst());
    
  }

  // MANUAL CHANGE
  void print(std::ostream &output_stream) const {
    output_stream << "Fo_Scalarss_Spherfuncs:" << std::endl
      << " starter: " << starter()  << std::endl
      << " maxmoms: " << maxmoms()  << std::endl
      << "maxgeoms: " << maxgeoms()  << std::endl
      << "   nmoms: " << nmoms()  << std::endl
      << "  ngeoms: " << ngeoms()  << std::endl
      << "     df1: " << std::endl << df1()  << std::endl
      << "     df2: " << std::endl << df2()  << std::endl
      << " cosscat: " << std::endl << cosscat()  << std::endl
      << " ss_pleg: " << std::endl << ss_pleg()  << std::endl;
  }
  // MANUAL CHANGE

private:
  bool starter_;
  int maxmoms_;
  int maxgeoms_;
  int nmoms_;
  int ngeoms_;
  blitz::Array<double, 1> df1_;
  blitz::Array<double, 1> df2_;
  blitz::Array<double, 1> cosscat_;
  blitz::Array<double, 2> ss_pleg_;
  // MANUAL CHANGE
  Fo_Scalarss_Spherfuncs() {}
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const;
  template<class Archive>
  void load(Archive & ar, const unsigned int version);
  // MANUAL CHANGE
};

//-----------------------------------------------------------------------
// Links to module: "fo_thermal_rtcalcs_i_m" in file: "FO_Thermal_RTCalcs_I.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_thermal_rtcalcs_i_m_dte_integral_i_dn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxfinelayers_in, const int* max_user_levels_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_planpar_in, const bool* do_regular_ps_in, const bool* do_enhanced_ps_in, const bool* donadir_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const double* bb_input_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* mu1_in, const int* ncrit_in, const double* radcrit_in, const double* cotcrit_in, const double* raycon_in, const double* radii_in, const double* cota_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* intensity_dta_dn_in, const double* cumsource_dn_in, const double* tcom1_in);
  void fo_thermal_rtcalcs_i_m_dte_integral_i_up_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxfinelayers_in, const int* max_user_levels_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_planpar_in, const bool* do_regular_ps_in, const bool* do_enhanced_ps_in, const bool* donadir_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const double* bb_input_in, const double* surfbb_in, const double* user_emissivity_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* mu1_in, const int* ncrit_in, const double* raycon_in, const double* cota_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* intensity_dta_up_in, const double* intensity_dts_in, const double* cumsource_up_in, const double* tcom1_in);
  void fo_thermal_rtcalcs_i_m_dte_integral_i_updn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxfinelayers_in, const int* max_user_levels_in, const bool* do_upwelling_in, const bool* do_dnwelling_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_planpar_in, const bool* do_regular_ps_in, const bool* do_enhanced_ps_in, const bool* donadir_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const double* bb_input_in, const double* surfbb_in, const double* user_emissivity_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* mu1_in, const int* ncrit_in, const double* radcrit_in, const double* cotcrit_in, const double* raycon_in, const double* radii_in, const double* cota_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* intensity_dta_up_in, const double* intensity_dts_in, const double* cumsource_up_in, const double* intensity_dta_dn_in, const double* cumsource_dn_in, const double* tcom1_in);
}

class Fo_Thermal_Rtcalcs_I : public virtual GenericObject {

public:
  Fo_Thermal_Rtcalcs_I(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfinelayers_in, const int& max_user_levels_in, const int& ngeoms_in, const int& nlayers_in, const int& n_user_levels_in) : maxgeoms_(maxgeoms_in), maxlayers_(maxlayers_in), maxfinelayers_(maxfinelayers_in), max_user_levels_(max_user_levels_in), ngeoms_(ngeoms_in), nlayers_(nlayers_in), n_user_levels_(n_user_levels_in) 
  { 
    do_thermset_ = false;
    do_deltam_scaling_ = false;
    do_planpar_ = false;
    do_regular_ps_ = false;
    do_enhanced_ps_ = false;
    donadir_.reference( blitz::Array<bool, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    donadir_ = false;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    user_levels_.reference( blitz::Array<int, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    user_levels_ = 0;
    bb_input_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    bb_input_ = 0;
    extinction_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    extinction_ = 0;
    deltaus_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    deltaus_ = 0;
    omega_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    omega_ = 0;
    truncfac_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    truncfac_ = 0;
    mu1_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu1_ = 0;
    ncrit_.reference( blitz::Array<int, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    ncrit_ = 0;
    radcrit_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    radcrit_ = 0;
    cotcrit_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    cotcrit_ = 0;
    raycon_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    raycon_ = 0;
    radii_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    radii_ = 0;
    cota_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cota_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfinelayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfinelayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    csqfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfinelayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    csqfine_ = 0;
    cotfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfinelayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    cotfine_ = 0;
    intensity_dta_dn_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dta_dn_ = 0;
    cumsource_dn_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cumsource_dn_ = 0;
    tcom1_.reference( blitz::Array<double, 2>(maxlayers_, 2, blitz::ColumnMajorArray<2>()) );
    tcom1_ = 0;
    surfbb_ = 0;
    user_emissivity_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    user_emissivity_ = 0;
    intensity_dta_up_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dta_up_ = 0;
    intensity_dts_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dts_ = 0;
    cumsource_up_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cumsource_up_ = 0;
    do_upwelling_ = false;
    do_dnwelling_ = false;
    // Initialize type pointers
    
  }

  virtual ~Fo_Thermal_Rtcalcs_I() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxfinelayers() const {
    return maxfinelayers_;
  }

  

  const int& max_user_levels() const {
    return max_user_levels_;
  }

  

  const bool& do_thermset() const {
    return do_thermset_;
  }

  void do_thermset(const bool& do_thermset_in) {
    do_thermset_ = do_thermset_in;
  }

  

  const bool& do_deltam_scaling() const {
    return do_deltam_scaling_;
  }

  void do_deltam_scaling(const bool& do_deltam_scaling_in) {
    do_deltam_scaling_ = do_deltam_scaling_in;
  }

  

  const bool& do_planpar() const {
    return do_planpar_;
  }

  void do_planpar(const bool& do_planpar_in) {
    do_planpar_ = do_planpar_in;
  }

  

  const bool& do_regular_ps() const {
    return do_regular_ps_;
  }

  void do_regular_ps(const bool& do_regular_ps_in) {
    do_regular_ps_ = do_regular_ps_in;
  }

  

  const bool& do_enhanced_ps() const {
    return do_enhanced_ps_;
  }

  void do_enhanced_ps(const bool& do_enhanced_ps_in) {
    do_enhanced_ps_ = do_enhanced_ps_in;
  }

  

  const blitz::Array<bool, 1>& donadir() const {
    return donadir_;
  }

  void donadir(const blitz::Array<bool, 1>& donadir_in) {
    donadir_ = donadir_in;
  }

  

  const int& ngeoms() const {
    return ngeoms_;
  }

  

  const int& nlayers() const {
    return nlayers_;
  }

  

  const blitz::Array<int, 2>& nfinedivs() const {
    return nfinedivs_;
  }

  void nfinedivs(const blitz::Array<int, 2>& nfinedivs_in) {
    nfinedivs_ = nfinedivs_in;
  }

  

  const int& n_user_levels() const {
    return n_user_levels_;
  }

  

  const blitz::Array<int, 1>& user_levels() const {
    return user_levels_;
  }

  void user_levels(const blitz::Array<int, 1>& user_levels_in) {
    user_levels_ = user_levels_in;
  }

  

  const blitz::Array<double, 1>& bb_input() const {
    return bb_input_;
  }

  void bb_input(const blitz::Array<double, 1>& bb_input_in) {
    bb_input_ = bb_input_in;
  }

  

  const blitz::Array<double, 1>& extinction() const {
    return extinction_;
  }

  void extinction(const blitz::Array<double, 1>& extinction_in) {
    extinction_ = extinction_in;
  }

  

  const blitz::Array<double, 1>& deltaus() const {
    return deltaus_;
  }

  void deltaus(const blitz::Array<double, 1>& deltaus_in) {
    deltaus_ = deltaus_in;
  }

  

  const blitz::Array<double, 1>& omega() const {
    return omega_;
  }

  void omega(const blitz::Array<double, 1>& omega_in) {
    omega_ = omega_in;
  }

  

  const blitz::Array<double, 1>& truncfac() const {
    return truncfac_;
  }

  void truncfac(const blitz::Array<double, 1>& truncfac_in) {
    truncfac_ = truncfac_in;
  }

  

  const blitz::Array<double, 1>& mu1() const {
    return mu1_;
  }

  void mu1(const blitz::Array<double, 1>& mu1_in) {
    mu1_ = mu1_in;
  }

  

  const blitz::Array<int, 1>& ncrit() const {
    return ncrit_;
  }

  void ncrit(const blitz::Array<int, 1>& ncrit_in) {
    ncrit_ = ncrit_in;
  }

  

  const blitz::Array<double, 1>& radcrit() const {
    return radcrit_;
  }

  void radcrit(const blitz::Array<double, 1>& radcrit_in) {
    radcrit_ = radcrit_in;
  }

  

  const blitz::Array<double, 1>& cotcrit() const {
    return cotcrit_;
  }

  void cotcrit(const blitz::Array<double, 1>& cotcrit_in) {
    cotcrit_ = cotcrit_in;
  }

  

  const blitz::Array<double, 1>& raycon() const {
    return raycon_;
  }

  void raycon(const blitz::Array<double, 1>& raycon_in) {
    raycon_ = raycon_in;
  }

  

  const blitz::Array<double, 1>& radii() const {
    return radii_;
  }

  void radii(const blitz::Array<double, 1>& radii_in) {
    radii_ = radii_in;
  }

  

  const blitz::Array<double, 2>& cota() const {
    return cota_;
  }

  void cota(const blitz::Array<double, 2>& cota_in) {
    cota_ = cota_in;
  }

  

  const blitz::Array<double, 3>& xfine() const {
    return xfine_;
  }

  void xfine(const blitz::Array<double, 3>& xfine_in) {
    xfine_ = xfine_in;
  }

  

  const blitz::Array<double, 3>& wfine() const {
    return wfine_;
  }

  void wfine(const blitz::Array<double, 3>& wfine_in) {
    wfine_ = wfine_in;
  }

  

  const blitz::Array<double, 3>& csqfine() const {
    return csqfine_;
  }

  void csqfine(const blitz::Array<double, 3>& csqfine_in) {
    csqfine_ = csqfine_in;
  }

  

  const blitz::Array<double, 3>& cotfine() const {
    return cotfine_;
  }

  void cotfine(const blitz::Array<double, 3>& cotfine_in) {
    cotfine_ = cotfine_in;
  }

  

  const blitz::Array<double, 2>& intensity_dta_dn() const {
    return intensity_dta_dn_;
  }

  

  const blitz::Array<double, 2>& cumsource_dn() const {
    return cumsource_dn_;
  }

  

  const blitz::Array<double, 2>& tcom1() const {
    return tcom1_;
  }

  void tcom1(const blitz::Array<double, 2>& tcom1_in) {
    tcom1_ = tcom1_in;
  }

  

  const double& surfbb() const {
    return surfbb_;
  }

  void surfbb(const double& surfbb_in) {
    surfbb_ = surfbb_in;
  }

  

  const blitz::Array<double, 1>& user_emissivity() const {
    return user_emissivity_;
  }

  void user_emissivity(const blitz::Array<double, 1>& user_emissivity_in) {
    user_emissivity_ = user_emissivity_in;
  }

  

  const blitz::Array<double, 2>& intensity_dta_up() const {
    return intensity_dta_up_;
  }

  

  const blitz::Array<double, 2>& intensity_dts() const {
    return intensity_dts_;
  }

  

  const blitz::Array<double, 2>& cumsource_up() const {
    return cumsource_up_;
  }

  

  const bool& do_upwelling() const {
    return do_upwelling_;
  }

  void do_upwelling(const bool& do_upwelling_in) {
    do_upwelling_ = do_upwelling_in;
  }

  

  const bool& do_dnwelling() const {
    return do_dnwelling_;
  }

  void do_dnwelling(const bool& do_dnwelling_in) {
    do_dnwelling_ = do_dnwelling_in;
  }

  

  
  void dte_integral_i_dn() {
    
    
    fo_thermal_rtcalcs_i_m_dte_integral_i_dn_wrap(&maxgeoms_, &maxlayers_, &maxfinelayers_, &max_user_levels_, &do_thermset_, &do_deltam_scaling_, &do_planpar_, &do_regular_ps_, &do_enhanced_ps_, donadir_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), bb_input_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), mu1_.dataFirst(), ncrit_.dataFirst(), radcrit_.dataFirst(), cotcrit_.dataFirst(), raycon_.dataFirst(), radii_.dataFirst(), cota_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), csqfine_.dataFirst(), cotfine_.dataFirst(), intensity_dta_dn_.dataFirst(), cumsource_dn_.dataFirst(), tcom1_.dataFirst());
    
  }
void dte_integral_i_up() {
    
    
    fo_thermal_rtcalcs_i_m_dte_integral_i_up_wrap(&maxgeoms_, &maxlayers_, &maxfinelayers_, &max_user_levels_, &do_thermset_, &do_deltam_scaling_, &do_planpar_, &do_regular_ps_, &do_enhanced_ps_, donadir_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), bb_input_.dataFirst(), &surfbb_, user_emissivity_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), mu1_.dataFirst(), ncrit_.dataFirst(), raycon_.dataFirst(), cota_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), csqfine_.dataFirst(), cotfine_.dataFirst(), intensity_dta_up_.dataFirst(), intensity_dts_.dataFirst(), cumsource_up_.dataFirst(), tcom1_.dataFirst());
    
  }
void dte_integral_i_updn() {
    
    
    fo_thermal_rtcalcs_i_m_dte_integral_i_updn_wrap(&maxgeoms_, &maxlayers_, &maxfinelayers_, &max_user_levels_, &do_upwelling_, &do_dnwelling_, &do_thermset_, &do_deltam_scaling_, &do_planpar_, &do_regular_ps_, &do_enhanced_ps_, donadir_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), bb_input_.dataFirst(), &surfbb_, user_emissivity_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), mu1_.dataFirst(), ncrit_.dataFirst(), radcrit_.dataFirst(), cotcrit_.dataFirst(), raycon_.dataFirst(), radii_.dataFirst(), cota_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), csqfine_.dataFirst(), cotfine_.dataFirst(), intensity_dta_up_.dataFirst(), intensity_dts_.dataFirst(), cumsource_up_.dataFirst(), intensity_dta_dn_.dataFirst(), cumsource_dn_.dataFirst(), tcom1_.dataFirst());
    
  }

  friend std::ostream& operator<<(std::ostream &output_stream, const Fo_Thermal_Rtcalcs_I &obj) {
    output_stream << "Fo_Thermal_Rtcalcs_I:" << std::endl
      << "         maxgeoms: " << obj.maxgeoms()  << std::endl
      << "        maxlayers: " << obj.maxlayers()  << std::endl
      << "    maxfinelayers: " << obj.maxfinelayers()  << std::endl
      << "  max_user_levels: " << obj.max_user_levels()  << std::endl
      << "      do_thermset: " << obj.do_thermset()  << std::endl
      << "do_deltam_scaling: " << obj.do_deltam_scaling()  << std::endl
      << "       do_planpar: " << obj.do_planpar()  << std::endl
      << "    do_regular_ps: " << obj.do_regular_ps()  << std::endl
      << "   do_enhanced_ps: " << obj.do_enhanced_ps()  << std::endl
      << "          donadir: " << std::endl << obj.donadir()  << std::endl
      << "           ngeoms: " << obj.ngeoms()  << std::endl
      << "          nlayers: " << obj.nlayers()  << std::endl
      << "        nfinedivs: " << std::endl << obj.nfinedivs()  << std::endl
      << "    n_user_levels: " << obj.n_user_levels()  << std::endl
      << "      user_levels: " << std::endl << obj.user_levels()  << std::endl
      << "         bb_input: " << std::endl << obj.bb_input()  << std::endl
      << "       extinction: " << std::endl << obj.extinction()  << std::endl
      << "          deltaus: " << std::endl << obj.deltaus()  << std::endl
      << "            omega: " << std::endl << obj.omega()  << std::endl
      << "         truncfac: " << std::endl << obj.truncfac()  << std::endl
      << "              mu1: " << std::endl << obj.mu1()  << std::endl
      << "            ncrit: " << std::endl << obj.ncrit()  << std::endl
      << "          radcrit: " << std::endl << obj.radcrit()  << std::endl
      << "          cotcrit: " << std::endl << obj.cotcrit()  << std::endl
      << "           raycon: " << std::endl << obj.raycon()  << std::endl
      << "            radii: " << std::endl << obj.radii()  << std::endl
      << "             cota: " << std::endl << obj.cota()  << std::endl
      << "            xfine: " << std::endl << obj.xfine()  << std::endl
      << "            wfine: " << std::endl << obj.wfine()  << std::endl
      << "          csqfine: " << std::endl << obj.csqfine()  << std::endl
      << "          cotfine: " << std::endl << obj.cotfine()  << std::endl
      << " intensity_dta_dn: " << std::endl << obj.intensity_dta_dn()  << std::endl
      << "     cumsource_dn: " << std::endl << obj.cumsource_dn()  << std::endl
      << "            tcom1: " << std::endl << obj.tcom1()  << std::endl
      << "           surfbb: " << obj.surfbb()  << std::endl
      << "  user_emissivity: " << std::endl << obj.user_emissivity()  << std::endl
      << " intensity_dta_up: " << std::endl << obj.intensity_dta_up()  << std::endl
      << "    intensity_dts: " << std::endl << obj.intensity_dts()  << std::endl
      << "     cumsource_up: " << std::endl << obj.cumsource_up()  << std::endl
      << "     do_upwelling: " << obj.do_upwelling()  << std::endl
      << "     do_dnwelling: " << obj.do_dnwelling()  << std::endl;
    return output_stream;

  }

private:
  int maxgeoms_;
  int maxlayers_;
  int maxfinelayers_;
  int max_user_levels_;
  bool do_thermset_;
  bool do_deltam_scaling_;
  bool do_planpar_;
  bool do_regular_ps_;
  bool do_enhanced_ps_;
  blitz::Array<bool, 1> donadir_;
  int ngeoms_;
  int nlayers_;
  blitz::Array<int, 2> nfinedivs_;
  int n_user_levels_;
  blitz::Array<int, 1> user_levels_;
  blitz::Array<double, 1> bb_input_;
  blitz::Array<double, 1> extinction_;
  blitz::Array<double, 1> deltaus_;
  blitz::Array<double, 1> omega_;
  blitz::Array<double, 1> truncfac_;
  blitz::Array<double, 1> mu1_;
  blitz::Array<int, 1> ncrit_;
  blitz::Array<double, 1> radcrit_;
  blitz::Array<double, 1> cotcrit_;
  blitz::Array<double, 1> raycon_;
  blitz::Array<double, 1> radii_;
  blitz::Array<double, 2> cota_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> csqfine_;
  blitz::Array<double, 3> cotfine_;
  blitz::Array<double, 2> intensity_dta_dn_;
  blitz::Array<double, 2> cumsource_dn_;
  blitz::Array<double, 2> tcom1_;
  double surfbb_;
  blitz::Array<double, 1> user_emissivity_;
  blitz::Array<double, 2> intensity_dta_up_;
  blitz::Array<double, 2> intensity_dts_;
  blitz::Array<double, 2> cumsource_up_;
  bool do_upwelling_;
  bool do_dnwelling_;
};

//-----------------------------------------------------------------------
// Links to module: "fo_thermal_rtcalcs_ilpsb_m" in file: "FO_Thermal_RTCalcs_ILPSB.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_thermal_rtcalcs_ilpsb_m_dte_integral_ilpsb_dn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxfinelayers_in, const int* max_user_levels_in, const int* max_atmoswfs_in, const bool* do_abbwf_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_planpar_in, const bool* do_regular_ps_in, const bool* do_enhanced_ps_in, const bool* donadir_in, const bool* do_profilewfs_in, const bool* lvaryflags_in, const int* lvarynums_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const double* bb_input_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* l_extinction_in, const double* l_deltaus_in, const double* l_omega_in, const double* l_truncfac_in, const double* mu1_in, const int* ncrit_in, const double* radcrit_in, const double* cotcrit_in, const double* raycon_in, const double* cota_in, const double* radii_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* intensity_dta_dn_in, const double* lp_jacobians_dta_dn_in, const double* lab_jacobians_dta_dn_in, const double* tcom_in, const double* l_tcom_in, const double* lb_tcom1_in, const double* lb_tcom2_in);
  void fo_thermal_rtcalcs_ilpsb_m_dte_integral_ilpsb_up_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxfinelayers_in, const int* max_user_levels_in, const int* max_atmoswfs_in, const int* max_surfacewfs_in, const bool* do_abbwf_in, const bool* do_sbbwf_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_planpar_in, const bool* do_regular_ps_in, const bool* do_enhanced_ps_in, const bool* donadir_in, const bool* do_profilewfs_in, const bool* do_surfacewfs_in, const bool* lvaryflags_in, const int* lvarynums_in, const int* n_surfacewfs_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const double* bb_input_in, const double* surfbb_in, const double* user_emissivity_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* ls_user_emissivity_in, const double* l_extinction_in, const double* l_deltaus_in, const double* l_omega_in, const double* l_truncfac_in, const double* mu1_in, const int* ncrit_in, const double* raycon_in, const double* cota_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* intensity_dta_up_in, const double* intensity_dts_in, const double* lp_jacobians_dta_up_in, const double* lp_jacobians_dts_up_in, const double* ls_jacobians_dts_in, const double* lab_jacobians_dta_up_in, const double* lsb_jacobians_dts_in, const double* tcom_in, const double* l_tcom_in, const double* lb_tcom1_in, const double* lb_tcom2_in);
  void fo_thermal_rtcalcs_ilpsb_m_dte_integral_ilpsb_updn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxfinelayers_in, const int* max_user_levels_in, const int* max_atmoswfs_in, const int* max_surfacewfs_in, const bool* do_abbwf_in, const bool* do_sbbwf_in, const bool* do_upwelling_in, const bool* do_dnwelling_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_planpar_in, const bool* do_regular_ps_in, const bool* do_enhanced_ps_in, const bool* donadir_in, const bool* do_profilewfs_in, const bool* do_surfacewfs_in, const bool* lvaryflags_in, const int* lvarynums_in, const int* n_surfacewfs_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const double* bb_input_in, const double* surfbb_in, const double* user_emissivity_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* ls_user_emissivity_in, const double* l_extinction_in, const double* l_deltaus_in, const double* l_omega_in, const double* l_truncfac_in, const double* mu1_in, const int* ncrit_in, const double* radcrit_in, const double* cotcrit_in, const double* raycon_in, const double* cota_in, const double* radii_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* intensity_dta_up_in, const double* intensity_dts_in, const double* intensity_dta_dn_in, const double* lp_jacobians_dta_up_in, const double* lp_jacobians_dts_up_in, const double* lab_jacobians_dta_up_in, const double* ls_jacobians_dts_in, const double* lsb_jacobians_dts_in, const double* lp_jacobians_dta_dn_in, const double* lab_jacobians_dta_dn_in, const double* tcom_in, const double* l_tcom_in, const double* lb_tcom1_in, const double* lb_tcom2_in);
}

class Fo_Thermal_Rtcalcs_Ilpsb : public virtual GenericObject {

public:
  Fo_Thermal_Rtcalcs_Ilpsb(const int& maxgeoms_in, const int& maxlayers_in, const int& maxfinelayers_in, const int& max_user_levels_in, const int& max_atmoswfs_in, const int& ngeoms_in, const int& nlayers_in, const int& n_user_levels_in, const int& max_surfacewfs_in) : maxgeoms_(maxgeoms_in), maxlayers_(maxlayers_in), maxfinelayers_(maxfinelayers_in), max_user_levels_(max_user_levels_in), max_atmoswfs_(max_atmoswfs_in), ngeoms_(ngeoms_in), nlayers_(nlayers_in), n_user_levels_(n_user_levels_in), max_surfacewfs_(max_surfacewfs_in) 
  { 
    do_abbwf_ = false;
    do_thermset_ = false;
    do_deltam_scaling_ = false;
    do_planpar_ = false;
    do_regular_ps_ = false;
    do_enhanced_ps_ = false;
    donadir_.reference( blitz::Array<bool, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    donadir_ = false;
    do_profilewfs_ = false;
    lvaryflags_.reference( blitz::Array<bool, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    lvaryflags_ = false;
    lvarynums_.reference( blitz::Array<int, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    lvarynums_ = 0;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    user_levels_.reference( blitz::Array<int, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    user_levels_ = 0;
    bb_input_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    bb_input_ = 0;
    extinction_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    extinction_ = 0;
    deltaus_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    deltaus_ = 0;
    omega_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    omega_ = 0;
    truncfac_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    truncfac_ = 0;
    l_extinction_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_extinction_ = 0;
    l_deltaus_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_deltaus_ = 0;
    l_omega_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_omega_ = 0;
    l_truncfac_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_truncfac_ = 0;
    mu1_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu1_ = 0;
    ncrit_.reference( blitz::Array<int, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    ncrit_ = 0;
    radcrit_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    radcrit_ = 0;
    cotcrit_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    cotcrit_ = 0;
    raycon_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    raycon_ = 0;
    cota_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cota_ = 0;
    radii_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    radii_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfinelayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfinelayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    csqfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfinelayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    csqfine_ = 0;
    cotfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfinelayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    cotfine_ = 0;
    intensity_dta_dn_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dta_dn_ = 0;
    lp_jacobians_dta_dn_.reference( blitz::Array<double, 4>(max_user_levels_, maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    lp_jacobians_dta_dn_ = 0;
    lab_jacobians_dta_dn_.reference( blitz::Array<double, 3>(max_user_levels_, maxgeoms_, maxlayers_-0+1, blitz::ColumnMajorArray<3>()) );
    lab_jacobians_dta_dn_ = 0;
    tcom_.reference( blitz::Array<double, 2>(maxlayers_, 2, blitz::ColumnMajorArray<2>()) );
    tcom_ = 0;
    l_tcom_.reference( blitz::Array<double, 3>(maxlayers_, 2, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    l_tcom_ = 0;
    lb_tcom1_.reference( blitz::Array<double, 2>(maxlayers_, 2, blitz::ColumnMajorArray<2>()) );
    lb_tcom1_ = 0;
    lb_tcom2_.reference( blitz::Array<double, 2>(maxlayers_, 2, blitz::ColumnMajorArray<2>()) );
    lb_tcom2_ = 0;
    do_sbbwf_ = false;
    do_surfacewfs_ = false;
    n_surfacewfs_ = 0;
    surfbb_ = 0;
    user_emissivity_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    user_emissivity_ = 0;
    ls_user_emissivity_.reference( blitz::Array<double, 2>(maxgeoms_, max_surfacewfs_, blitz::ColumnMajorArray<2>()) );
    ls_user_emissivity_ = 0;
    intensity_dta_up_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dta_up_ = 0;
    intensity_dts_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dts_ = 0;
    lp_jacobians_dta_up_.reference( blitz::Array<double, 4>(max_user_levels_, maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    lp_jacobians_dta_up_ = 0;
    lp_jacobians_dts_up_.reference( blitz::Array<double, 4>(max_user_levels_, maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    lp_jacobians_dts_up_ = 0;
    ls_jacobians_dts_.reference( blitz::Array<double, 3>(max_user_levels_, maxgeoms_, max_surfacewfs_, blitz::ColumnMajorArray<3>()) );
    ls_jacobians_dts_ = 0;
    lab_jacobians_dta_up_.reference( blitz::Array<double, 3>(max_user_levels_, maxgeoms_, maxlayers_-0+1, blitz::ColumnMajorArray<3>()) );
    lab_jacobians_dta_up_ = 0;
    lsb_jacobians_dts_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    lsb_jacobians_dts_ = 0;
    do_upwelling_ = false;
    do_dnwelling_ = false;
    // Initialize type pointers
    
  }

  virtual ~Fo_Thermal_Rtcalcs_Ilpsb() = default;

  std::string print_to_string() const
  {
      std::ostringstream output;
      output << *this;
      return output.str();
  }

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxfinelayers() const {
    return maxfinelayers_;
  }

  

  const int& max_user_levels() const {
    return max_user_levels_;
  }

  

  const int& max_atmoswfs() const {
    return max_atmoswfs_;
  }

  

  const bool& do_abbwf() const {
    return do_abbwf_;
  }

  void do_abbwf(const bool& do_abbwf_in) {
    do_abbwf_ = do_abbwf_in;
  }

  

  const bool& do_thermset() const {
    return do_thermset_;
  }

  void do_thermset(const bool& do_thermset_in) {
    do_thermset_ = do_thermset_in;
  }

  

  const bool& do_deltam_scaling() const {
    return do_deltam_scaling_;
  }

  void do_deltam_scaling(const bool& do_deltam_scaling_in) {
    do_deltam_scaling_ = do_deltam_scaling_in;
  }

  

  const bool& do_planpar() const {
    return do_planpar_;
  }

  void do_planpar(const bool& do_planpar_in) {
    do_planpar_ = do_planpar_in;
  }

  

  const bool& do_regular_ps() const {
    return do_regular_ps_;
  }

  void do_regular_ps(const bool& do_regular_ps_in) {
    do_regular_ps_ = do_regular_ps_in;
  }

  

  const bool& do_enhanced_ps() const {
    return do_enhanced_ps_;
  }

  void do_enhanced_ps(const bool& do_enhanced_ps_in) {
    do_enhanced_ps_ = do_enhanced_ps_in;
  }

  

  const blitz::Array<bool, 1>& donadir() const {
    return donadir_;
  }

  void donadir(const blitz::Array<bool, 1>& donadir_in) {
    donadir_ = donadir_in;
  }

  

  const bool& do_profilewfs() const {
    return do_profilewfs_;
  }

  void do_profilewfs(const bool& do_profilewfs_in) {
    do_profilewfs_ = do_profilewfs_in;
  }

  

  const blitz::Array<bool, 1>& lvaryflags() const {
    return lvaryflags_;
  }

  void lvaryflags(const blitz::Array<bool, 1>& lvaryflags_in) {
    lvaryflags_ = lvaryflags_in;
  }

  

  const blitz::Array<int, 1>& lvarynums() const {
    return lvarynums_;
  }

  void lvarynums(const blitz::Array<int, 1>& lvarynums_in) {
    lvarynums_ = lvarynums_in;
  }

  

  const int& ngeoms() const {
    return ngeoms_;
  }

  

  const int& nlayers() const {
    return nlayers_;
  }

  

  const blitz::Array<int, 2>& nfinedivs() const {
    return nfinedivs_;
  }

  void nfinedivs(const blitz::Array<int, 2>& nfinedivs_in) {
    nfinedivs_ = nfinedivs_in;
  }

  

  const int& n_user_levels() const {
    return n_user_levels_;
  }

  

  const blitz::Array<int, 1>& user_levels() const {
    return user_levels_;
  }

  void user_levels(const blitz::Array<int, 1>& user_levels_in) {
    user_levels_ = user_levels_in;
  }

  

  const blitz::Array<double, 1>& bb_input() const {
    return bb_input_;
  }

  void bb_input(const blitz::Array<double, 1>& bb_input_in) {
    bb_input_ = bb_input_in;
  }

  

  const blitz::Array<double, 1>& extinction() const {
    return extinction_;
  }

  void extinction(const blitz::Array<double, 1>& extinction_in) {
    extinction_ = extinction_in;
  }

  

  const blitz::Array<double, 1>& deltaus() const {
    return deltaus_;
  }

  void deltaus(const blitz::Array<double, 1>& deltaus_in) {
    deltaus_ = deltaus_in;
  }

  

  const blitz::Array<double, 1>& omega() const {
    return omega_;
  }

  void omega(const blitz::Array<double, 1>& omega_in) {
    omega_ = omega_in;
  }

  

  const blitz::Array<double, 1>& truncfac() const {
    return truncfac_;
  }

  void truncfac(const blitz::Array<double, 1>& truncfac_in) {
    truncfac_ = truncfac_in;
  }

  

  const blitz::Array<double, 2>& l_extinction() const {
    return l_extinction_;
  }

  void l_extinction(const blitz::Array<double, 2>& l_extinction_in) {
    l_extinction_ = l_extinction_in;
  }

  

  const blitz::Array<double, 2>& l_deltaus() const {
    return l_deltaus_;
  }

  void l_deltaus(const blitz::Array<double, 2>& l_deltaus_in) {
    l_deltaus_ = l_deltaus_in;
  }

  

  const blitz::Array<double, 2>& l_omega() const {
    return l_omega_;
  }

  void l_omega(const blitz::Array<double, 2>& l_omega_in) {
    l_omega_ = l_omega_in;
  }

  

  const blitz::Array<double, 2>& l_truncfac() const {
    return l_truncfac_;
  }

  void l_truncfac(const blitz::Array<double, 2>& l_truncfac_in) {
    l_truncfac_ = l_truncfac_in;
  }

  

  const blitz::Array<double, 1>& mu1() const {
    return mu1_;
  }

  void mu1(const blitz::Array<double, 1>& mu1_in) {
    mu1_ = mu1_in;
  }

  

  const blitz::Array<int, 1>& ncrit() const {
    return ncrit_;
  }

  void ncrit(const blitz::Array<int, 1>& ncrit_in) {
    ncrit_ = ncrit_in;
  }

  

  const blitz::Array<double, 1>& radcrit() const {
    return radcrit_;
  }

  void radcrit(const blitz::Array<double, 1>& radcrit_in) {
    radcrit_ = radcrit_in;
  }

  

  const blitz::Array<double, 1>& cotcrit() const {
    return cotcrit_;
  }

  void cotcrit(const blitz::Array<double, 1>& cotcrit_in) {
    cotcrit_ = cotcrit_in;
  }

  

  const blitz::Array<double, 1>& raycon() const {
    return raycon_;
  }

  void raycon(const blitz::Array<double, 1>& raycon_in) {
    raycon_ = raycon_in;
  }

  

  const blitz::Array<double, 2>& cota() const {
    return cota_;
  }

  void cota(const blitz::Array<double, 2>& cota_in) {
    cota_ = cota_in;
  }

  

  const blitz::Array<double, 1>& radii() const {
    return radii_;
  }

  void radii(const blitz::Array<double, 1>& radii_in) {
    radii_ = radii_in;
  }

  

  const blitz::Array<double, 3>& xfine() const {
    return xfine_;
  }

  void xfine(const blitz::Array<double, 3>& xfine_in) {
    xfine_ = xfine_in;
  }

  

  const blitz::Array<double, 3>& wfine() const {
    return wfine_;
  }

  void wfine(const blitz::Array<double, 3>& wfine_in) {
    wfine_ = wfine_in;
  }

  

  const blitz::Array<double, 3>& csqfine() const {
    return csqfine_;
  }

  void csqfine(const blitz::Array<double, 3>& csqfine_in) {
    csqfine_ = csqfine_in;
  }

  

  const blitz::Array<double, 3>& cotfine() const {
    return cotfine_;
  }

  void cotfine(const blitz::Array<double, 3>& cotfine_in) {
    cotfine_ = cotfine_in;
  }

  

  const blitz::Array<double, 2>& intensity_dta_dn() const {
    return intensity_dta_dn_;
  }

  

  const blitz::Array<double, 4>& lp_jacobians_dta_dn() const {
    return lp_jacobians_dta_dn_;
  }

  

  const blitz::Array<double, 3>& lab_jacobians_dta_dn() const {
    return lab_jacobians_dta_dn_;
  }

  

  const blitz::Array<double, 2>& tcom() const {
    return tcom_;
  }

  void tcom(const blitz::Array<double, 2>& tcom_in) {
    tcom_ = tcom_in;
  }

  

  const blitz::Array<double, 3>& l_tcom() const {
    return l_tcom_;
  }

  void l_tcom(const blitz::Array<double, 3>& l_tcom_in) {
    l_tcom_ = l_tcom_in;
  }

  

  const blitz::Array<double, 2>& lb_tcom1() const {
    return lb_tcom1_;
  }

  void lb_tcom1(const blitz::Array<double, 2>& lb_tcom1_in) {
    lb_tcom1_ = lb_tcom1_in;
  }

  

  const blitz::Array<double, 2>& lb_tcom2() const {
    return lb_tcom2_;
  }

  void lb_tcom2(const blitz::Array<double, 2>& lb_tcom2_in) {
    lb_tcom2_ = lb_tcom2_in;
  }

  

  const int& max_surfacewfs() const {
    return max_surfacewfs_;
  }

  

  const bool& do_sbbwf() const {
    return do_sbbwf_;
  }

  void do_sbbwf(const bool& do_sbbwf_in) {
    do_sbbwf_ = do_sbbwf_in;
  }

  

  const bool& do_surfacewfs() const {
    return do_surfacewfs_;
  }

  void do_surfacewfs(const bool& do_surfacewfs_in) {
    do_surfacewfs_ = do_surfacewfs_in;
  }

  

  const int& n_surfacewfs() const {
    return n_surfacewfs_;
  }

  void n_surfacewfs(const int& n_surfacewfs_in) {
    n_surfacewfs_ = n_surfacewfs_in;
  }

  

  const double& surfbb() const {
    return surfbb_;
  }

  void surfbb(const double& surfbb_in) {
    surfbb_ = surfbb_in;
  }

  

  const blitz::Array<double, 1>& user_emissivity() const {
    return user_emissivity_;
  }

  void user_emissivity(const blitz::Array<double, 1>& user_emissivity_in) {
    user_emissivity_ = user_emissivity_in;
  }

  

  const blitz::Array<double, 2>& ls_user_emissivity() const {
    return ls_user_emissivity_;
  }

  void ls_user_emissivity(const blitz::Array<double, 2>& ls_user_emissivity_in) {
    ls_user_emissivity_ = ls_user_emissivity_in;
  }

  

  const blitz::Array<double, 2>& intensity_dta_up() const {
    return intensity_dta_up_;
  }

  

  const blitz::Array<double, 2>& intensity_dts() const {
    return intensity_dts_;
  }

  

  const blitz::Array<double, 4>& lp_jacobians_dta_up() const {
    return lp_jacobians_dta_up_;
  }

  

  const blitz::Array<double, 4>& lp_jacobians_dts_up() const {
    return lp_jacobians_dts_up_;
  }

  

  const blitz::Array<double, 3>& ls_jacobians_dts() const {
    return ls_jacobians_dts_;
  }

  

  const blitz::Array<double, 3>& lab_jacobians_dta_up() const {
    return lab_jacobians_dta_up_;
  }

  

  const blitz::Array<double, 2>& lsb_jacobians_dts() const {
    return lsb_jacobians_dts_;
  }

  

  const bool& do_upwelling() const {
    return do_upwelling_;
  }

  void do_upwelling(const bool& do_upwelling_in) {
    do_upwelling_ = do_upwelling_in;
  }

  

  const bool& do_dnwelling() const {
    return do_dnwelling_;
  }

  void do_dnwelling(const bool& do_dnwelling_in) {
    do_dnwelling_ = do_dnwelling_in;
  }

  

  
  void dte_integral_ilpsb_dn() {
    
    
    fo_thermal_rtcalcs_ilpsb_m_dte_integral_ilpsb_dn_wrap(&maxgeoms_, &maxlayers_, &maxfinelayers_, &max_user_levels_, &max_atmoswfs_, &do_abbwf_, &do_thermset_, &do_deltam_scaling_, &do_planpar_, &do_regular_ps_, &do_enhanced_ps_, donadir_.dataFirst(), &do_profilewfs_, lvaryflags_.dataFirst(), lvarynums_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), bb_input_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), l_extinction_.dataFirst(), l_deltaus_.dataFirst(), l_omega_.dataFirst(), l_truncfac_.dataFirst(), mu1_.dataFirst(), ncrit_.dataFirst(), radcrit_.dataFirst(), cotcrit_.dataFirst(), raycon_.dataFirst(), cota_.dataFirst(), radii_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), csqfine_.dataFirst(), cotfine_.dataFirst(), intensity_dta_dn_.dataFirst(), lp_jacobians_dta_dn_.dataFirst(), lab_jacobians_dta_dn_.dataFirst(), tcom_.dataFirst(), l_tcom_.dataFirst(), lb_tcom1_.dataFirst(), lb_tcom2_.dataFirst());
    
  }
void dte_integral_ilpsb_up() {
    
    
    fo_thermal_rtcalcs_ilpsb_m_dte_integral_ilpsb_up_wrap(&maxgeoms_, &maxlayers_, &maxfinelayers_, &max_user_levels_, &max_atmoswfs_, &max_surfacewfs_, &do_abbwf_, &do_sbbwf_, &do_thermset_, &do_deltam_scaling_, &do_planpar_, &do_regular_ps_, &do_enhanced_ps_, donadir_.dataFirst(), &do_profilewfs_, &do_surfacewfs_, lvaryflags_.dataFirst(), lvarynums_.dataFirst(), &n_surfacewfs_, &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), bb_input_.dataFirst(), &surfbb_, user_emissivity_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), ls_user_emissivity_.dataFirst(), l_extinction_.dataFirst(), l_deltaus_.dataFirst(), l_omega_.dataFirst(), l_truncfac_.dataFirst(), mu1_.dataFirst(), ncrit_.dataFirst(), raycon_.dataFirst(), cota_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), csqfine_.dataFirst(), cotfine_.dataFirst(), intensity_dta_up_.dataFirst(), intensity_dts_.dataFirst(), lp_jacobians_dta_up_.dataFirst(), lp_jacobians_dts_up_.dataFirst(), ls_jacobians_dts_.dataFirst(), lab_jacobians_dta_up_.dataFirst(), lsb_jacobians_dts_.dataFirst(), tcom_.dataFirst(), l_tcom_.dataFirst(), lb_tcom1_.dataFirst(), lb_tcom2_.dataFirst());
    
  }
void dte_integral_ilpsb_updn() {
    
    
    fo_thermal_rtcalcs_ilpsb_m_dte_integral_ilpsb_updn_wrap(&maxgeoms_, &maxlayers_, &maxfinelayers_, &max_user_levels_, &max_atmoswfs_, &max_surfacewfs_, &do_abbwf_, &do_sbbwf_, &do_upwelling_, &do_dnwelling_, &do_thermset_, &do_deltam_scaling_, &do_planpar_, &do_regular_ps_, &do_enhanced_ps_, donadir_.dataFirst(), &do_profilewfs_, &do_surfacewfs_, lvaryflags_.dataFirst(), lvarynums_.dataFirst(), &n_surfacewfs_, &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), bb_input_.dataFirst(), &surfbb_, user_emissivity_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), ls_user_emissivity_.dataFirst(), l_extinction_.dataFirst(), l_deltaus_.dataFirst(), l_omega_.dataFirst(), l_truncfac_.dataFirst(), mu1_.dataFirst(), ncrit_.dataFirst(), radcrit_.dataFirst(), cotcrit_.dataFirst(), raycon_.dataFirst(), cota_.dataFirst(), radii_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), csqfine_.dataFirst(), cotfine_.dataFirst(), intensity_dta_up_.dataFirst(), intensity_dts_.dataFirst(), intensity_dta_dn_.dataFirst(), lp_jacobians_dta_up_.dataFirst(), lp_jacobians_dts_up_.dataFirst(), lab_jacobians_dta_up_.dataFirst(), ls_jacobians_dts_.dataFirst(), lsb_jacobians_dts_.dataFirst(), lp_jacobians_dta_dn_.dataFirst(), lab_jacobians_dta_dn_.dataFirst(), tcom_.dataFirst(), l_tcom_.dataFirst(), lb_tcom1_.dataFirst(), lb_tcom2_.dataFirst());
    
  }

  friend std::ostream& operator<<(std::ostream &output_stream, const Fo_Thermal_Rtcalcs_Ilpsb &obj) {
    output_stream << "Fo_Thermal_Rtcalcs_Ilpsb:" << std::endl
      << "            maxgeoms: " << obj.maxgeoms()  << std::endl
      << "           maxlayers: " << obj.maxlayers()  << std::endl
      << "       maxfinelayers: " << obj.maxfinelayers()  << std::endl
      << "     max_user_levels: " << obj.max_user_levels()  << std::endl
      << "        max_atmoswfs: " << obj.max_atmoswfs()  << std::endl
      << "            do_abbwf: " << obj.do_abbwf()  << std::endl
      << "         do_thermset: " << obj.do_thermset()  << std::endl
      << "   do_deltam_scaling: " << obj.do_deltam_scaling()  << std::endl
      << "          do_planpar: " << obj.do_planpar()  << std::endl
      << "       do_regular_ps: " << obj.do_regular_ps()  << std::endl
      << "      do_enhanced_ps: " << obj.do_enhanced_ps()  << std::endl
      << "             donadir: " << std::endl << obj.donadir()  << std::endl
      << "       do_profilewfs: " << obj.do_profilewfs()  << std::endl
      << "          lvaryflags: " << std::endl << obj.lvaryflags()  << std::endl
      << "           lvarynums: " << std::endl << obj.lvarynums()  << std::endl
      << "              ngeoms: " << obj.ngeoms()  << std::endl
      << "             nlayers: " << obj.nlayers()  << std::endl
      << "           nfinedivs: " << std::endl << obj.nfinedivs()  << std::endl
      << "       n_user_levels: " << obj.n_user_levels()  << std::endl
      << "         user_levels: " << std::endl << obj.user_levels()  << std::endl
      << "            bb_input: " << std::endl << obj.bb_input()  << std::endl
      << "          extinction: " << std::endl << obj.extinction()  << std::endl
      << "             deltaus: " << std::endl << obj.deltaus()  << std::endl
      << "               omega: " << std::endl << obj.omega()  << std::endl
      << "            truncfac: " << std::endl << obj.truncfac()  << std::endl
      << "        l_extinction: " << std::endl << obj.l_extinction()  << std::endl
      << "           l_deltaus: " << std::endl << obj.l_deltaus()  << std::endl
      << "             l_omega: " << std::endl << obj.l_omega()  << std::endl
      << "          l_truncfac: " << std::endl << obj.l_truncfac()  << std::endl
      << "                 mu1: " << std::endl << obj.mu1()  << std::endl
      << "               ncrit: " << std::endl << obj.ncrit()  << std::endl
      << "             radcrit: " << std::endl << obj.radcrit()  << std::endl
      << "             cotcrit: " << std::endl << obj.cotcrit()  << std::endl
      << "              raycon: " << std::endl << obj.raycon()  << std::endl
      << "                cota: " << std::endl << obj.cota()  << std::endl
      << "               radii: " << std::endl << obj.radii()  << std::endl
      << "               xfine: " << std::endl << obj.xfine()  << std::endl
      << "               wfine: " << std::endl << obj.wfine()  << std::endl
      << "             csqfine: " << std::endl << obj.csqfine()  << std::endl
      << "             cotfine: " << std::endl << obj.cotfine()  << std::endl
      << "    intensity_dta_dn: " << std::endl << obj.intensity_dta_dn()  << std::endl
      << " lp_jacobians_dta_dn: " << std::endl << obj.lp_jacobians_dta_dn()  << std::endl
      << "lab_jacobians_dta_dn: " << std::endl << obj.lab_jacobians_dta_dn()  << std::endl
      << "                tcom: " << std::endl << obj.tcom()  << std::endl
      << "              l_tcom: " << std::endl << obj.l_tcom()  << std::endl
      << "            lb_tcom1: " << std::endl << obj.lb_tcom1()  << std::endl
      << "            lb_tcom2: " << std::endl << obj.lb_tcom2()  << std::endl
      << "      max_surfacewfs: " << obj.max_surfacewfs()  << std::endl
      << "            do_sbbwf: " << obj.do_sbbwf()  << std::endl
      << "       do_surfacewfs: " << obj.do_surfacewfs()  << std::endl
      << "        n_surfacewfs: " << obj.n_surfacewfs()  << std::endl
      << "              surfbb: " << obj.surfbb()  << std::endl
      << "     user_emissivity: " << std::endl << obj.user_emissivity()  << std::endl
      << "  ls_user_emissivity: " << std::endl << obj.ls_user_emissivity()  << std::endl
      << "    intensity_dta_up: " << std::endl << obj.intensity_dta_up()  << std::endl
      << "       intensity_dts: " << std::endl << obj.intensity_dts()  << std::endl
      << " lp_jacobians_dta_up: " << std::endl << obj.lp_jacobians_dta_up()  << std::endl
      << " lp_jacobians_dts_up: " << std::endl << obj.lp_jacobians_dts_up()  << std::endl
      << "    ls_jacobians_dts: " << std::endl << obj.ls_jacobians_dts()  << std::endl
      << "lab_jacobians_dta_up: " << std::endl << obj.lab_jacobians_dta_up()  << std::endl
      << "   lsb_jacobians_dts: " << std::endl << obj.lsb_jacobians_dts()  << std::endl
      << "        do_upwelling: " << obj.do_upwelling()  << std::endl
      << "        do_dnwelling: " << obj.do_dnwelling()  << std::endl;
    return output_stream;

  }

private:
  int maxgeoms_;
  int maxlayers_;
  int maxfinelayers_;
  int max_user_levels_;
  int max_atmoswfs_;
  bool do_abbwf_;
  bool do_thermset_;
  bool do_deltam_scaling_;
  bool do_planpar_;
  bool do_regular_ps_;
  bool do_enhanced_ps_;
  blitz::Array<bool, 1> donadir_;
  bool do_profilewfs_;
  blitz::Array<bool, 1> lvaryflags_;
  blitz::Array<int, 1> lvarynums_;
  int ngeoms_;
  int nlayers_;
  blitz::Array<int, 2> nfinedivs_;
  int n_user_levels_;
  blitz::Array<int, 1> user_levels_;
  blitz::Array<double, 1> bb_input_;
  blitz::Array<double, 1> extinction_;
  blitz::Array<double, 1> deltaus_;
  blitz::Array<double, 1> omega_;
  blitz::Array<double, 1> truncfac_;
  blitz::Array<double, 2> l_extinction_;
  blitz::Array<double, 2> l_deltaus_;
  blitz::Array<double, 2> l_omega_;
  blitz::Array<double, 2> l_truncfac_;
  blitz::Array<double, 1> mu1_;
  blitz::Array<int, 1> ncrit_;
  blitz::Array<double, 1> radcrit_;
  blitz::Array<double, 1> cotcrit_;
  blitz::Array<double, 1> raycon_;
  blitz::Array<double, 2> cota_;
  blitz::Array<double, 1> radii_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> csqfine_;
  blitz::Array<double, 3> cotfine_;
  blitz::Array<double, 2> intensity_dta_dn_;
  blitz::Array<double, 4> lp_jacobians_dta_dn_;
  blitz::Array<double, 3> lab_jacobians_dta_dn_;
  blitz::Array<double, 2> tcom_;
  blitz::Array<double, 3> l_tcom_;
  blitz::Array<double, 2> lb_tcom1_;
  blitz::Array<double, 2> lb_tcom2_;
  int max_surfacewfs_;
  bool do_sbbwf_;
  bool do_surfacewfs_;
  int n_surfacewfs_;
  double surfbb_;
  blitz::Array<double, 1> user_emissivity_;
  blitz::Array<double, 2> ls_user_emissivity_;
  blitz::Array<double, 2> intensity_dta_up_;
  blitz::Array<double, 2> intensity_dts_;
  blitz::Array<double, 4> lp_jacobians_dta_up_;
  blitz::Array<double, 4> lp_jacobians_dts_up_;
  blitz::Array<double, 3> ls_jacobians_dts_;
  blitz::Array<double, 3> lab_jacobians_dta_up_;
  blitz::Array<double, 2> lsb_jacobians_dts_;
  bool do_upwelling_;
  bool do_dnwelling_;
};



}
  // MANUAL CHANGE
FP_EXPORT_KEY(Fo_Ssgeometry_Master)
FP_EXPORT_KEY(Fo_Scalarss_Spherfuncs)
FP_EXPORT_KEY(Fo_Scalarss_Rtcalcs_Ilps)
  // MANUAL CHANGE
#endif
