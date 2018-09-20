#ifndef FIRST_ORDER_INTERFACE_H
#define FIRST_ORDER_INTERFACE_H

#include <iostream>
#include <blitz/array.h>

#include "fp_exception.h"


/* This file was auto-generated */

namespace FullPhysics {

//-----------------------------------------------------------------------
// Links to module: "fo_dtgeometry_master_m" in file: "FO_DTgeometry_master.f90"
//-----------------------------------------------------------------------

extern "C" {

  void fo_dtgeometry_master_m_fo_dtgeometry_master_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxfine_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const int* ngeoms_in, const int* nlayers_in, const int* nfine_in, const double* dtr_in, const double* eradius_in, const double* heights_in, const double* alpha_boa_in, const bool* donadir_in, const bool* docrit_in, const double* acrit_in, const double* extinc_in, const double* raycon_in, const double* radii_in, const double* alpha_in, const double* cota_in, const int* nfinedivs_in, const double* xfine_in, const double* wfine_in, const double* csqfine_in, const double* cotfine_in, const double* alphafine_in, const double* radiifine_in, const int* ncrit_in, const double* radcrit_in, const double* cotcrit_in, const double* mu1_in, const bool* fail_in, const int* message_len, const char* message_in, const int* trace_len, const char* trace_in);
}

class Fo_Dtgeometry_Master {

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

class Fo_Ssgeometry_Master {

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

  friend std::ostream& operator<<(std::ostream &output_stream, const Fo_Ssgeometry_Master &obj) {
    output_stream << "Fo_Ssgeometry_Master:" << std::endl
      << "      maxgeoms: " << obj.maxgeoms()  << std::endl
      << "       maxszas: " << obj.maxszas()  << std::endl
      << "       maxvzas: " << obj.maxvzas()  << std::endl
      << "       maxazms: " << obj.maxazms()  << std::endl
      << "     maxlayers: " << obj.maxlayers()  << std::endl
      << "       maxfine: " << obj.maxfine()  << std::endl
      << "    do_obsgeom: " << obj.do_obsgeom()  << std::endl
      << "    do_chapman: " << obj.do_chapman()  << std::endl
      << "    do_planpar: " << obj.do_planpar()  << std::endl
      << "do_enhanced_ps: " << obj.do_enhanced_ps()  << std::endl
      << "        ngeoms: " << obj.ngeoms()  << std::endl
      << "         nszas: " << obj.nszas()  << std::endl
      << "         nvzas: " << obj.nvzas()  << std::endl
      << "         nazms: " << obj.nazms()  << std::endl
      << "       nlayers: " << obj.nlayers()  << std::endl
      << "         nfine: " << obj.nfine()  << std::endl
      << "           dtr: " << obj.dtr()  << std::endl
      << "           pie: " << obj.pie()  << std::endl
      << "         vsign: " << obj.vsign()  << std::endl
      << "       eradius: " << obj.eradius()  << std::endl
      << "       heights: " << std::endl << obj.heights()  << std::endl
      << "   obsgeom_boa: " << std::endl << obj.obsgeom_boa()  << std::endl
      << "     alpha_boa: " << std::endl << obj.alpha_boa()  << std::endl
      << "     theta_boa: " << std::endl << obj.theta_boa()  << std::endl
      << "       phi_boa: " << std::endl << obj.phi_boa()  << std::endl
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
      << "           mu0: " << std::endl << obj.mu0()  << std::endl
      << "           mu1: " << std::endl << obj.mu1()  << std::endl
      << "       cosscat: " << std::endl << obj.cosscat()  << std::endl
      << "      chapfacs: " << std::endl << obj.chapfacs()  << std::endl
      << "      sunpaths: " << std::endl << obj.sunpaths()  << std::endl
      << "     ntraverse: " << std::endl << obj.ntraverse()  << std::endl
      << "  sunpathsfine: " << std::endl << obj.sunpathsfine()  << std::endl
      << " ntraversefine: " << std::endl << obj.ntraversefine()  << std::endl
      << "          fail: " << obj.fail()  << std::endl
      << "       message: " << "\"" << obj.message() << "\"" << std::endl
      << "         trace: " << "\"" << obj.trace() << "\"" << std::endl;
    return output_stream;

  }

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
};



}
#endif