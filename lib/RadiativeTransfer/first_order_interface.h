#ifndef FIRST_ORDER_INTERFACE_H
#define FIRST_ORDER_INTERFACE_H

#include <iostream>
#include <blitz/array.h>

#include "fp_exception.h"


/* This file was auto-generated */

namespace FullPhysics {

//-----------------------------------------------------------------------
// Links to module: "fo_dtwpgeometry_master_m" in file: "FO_DTWPgeometry_master.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_dtwpgeometry_master_m_fo_dtwpgeometry_master_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const double* dtr_in, const double* eradius_in, const bool* do_upwelling_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const bool* do_partials_in, const int* ngeoms_in, const int* nlayers_in, const int* npartials_in, const int* nfine_in, const int* partial_layeridx_in, const double* heights_in, const double* alpha_boa_in, const double* partial_heights_in, const double* mu1_in, const double* radii_in, const double* losw_paths_in, const double* alpha_in, const double* sina_in, const double* cosa_in, const double* radii_p_in, const double* losp_paths_in, const double* alpha_p_in, const double* sina_p_in, const double* cosa_p_in, const int* nfinedivs_in, const double* xfine_in, const double* wfine_in, const double* radiifine_in, const double* alphafine_in, const double* sinfine_in, const double* cosfine_in, const int* nfinedivs_p_in, const double* xfine_p_in, const double* wfine_p_in, const double* radiifine_p_in, const double* alphafine_p_in, const double* sinfine_p_in, const double* cosfine_p_in, const bool* fail_in, const int* message_len, const char* message_in, const int* trace_len, const char* trace_in);
}

class Fo_Dtwpgeometry_Master : public Printable<Fo_Dtwpgeometry_Master> {

public:
  Fo_Dtwpgeometry_Master(const int& maxgeoms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& ngeoms_in, const int& nlayers_in, const int& npartials_in, const int& nfine_in) : maxgeoms_(maxgeoms_in), maxlayers_(maxlayers_in), maxpartials_(maxpartials_in), maxfine_(maxfine_in), ngeoms_(ngeoms_in), nlayers_(nlayers_in), npartials_(npartials_in), nfine_(nfine_in) 
  { 
    dtr_ = 0;
    eradius_ = 0;
    do_upwelling_ = false;
    do_planpar_ = false;
    do_enhanced_ps_ = false;
    do_partials_ = false;
    partial_layeridx_.reference( blitz::Array<int, 1>(maxpartials_, blitz::ColumnMajorArray<1>()) );
    partial_layeridx_ = 0;
    heights_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    heights_ = 0;
    alpha_boa_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    alpha_boa_ = 0;
    partial_heights_.reference( blitz::Array<double, 1>(maxpartials_, blitz::ColumnMajorArray<1>()) );
    partial_heights_ = 0;
    mu1_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu1_ = 0;
    radii_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    radii_ = 0;
    losw_paths_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losw_paths_ = 0;
    alpha_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    alpha_ = 0;
    sina_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    sina_ = 0;
    cosa_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cosa_ = 0;
    radii_p_.reference( blitz::Array<double, 1>(maxpartials_, blitz::ColumnMajorArray<1>()) );
    radii_p_ = 0;
    losp_paths_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losp_paths_ = 0;
    alpha_p_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    alpha_p_ = 0;
    sina_p_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    sina_p_ = 0;
    cosa_p_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cosa_p_ = 0;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    radiifine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    radiifine_ = 0;
    alphafine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    alphafine_ = 0;
    sinfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sinfine_ = 0;
    cosfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    cosfine_ = 0;
    nfinedivs_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_p_ = 0;
    xfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_p_ = 0;
    wfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_p_ = 0;
    radiifine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    radiifine_p_ = 0;
    alphafine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    alphafine_p_ = 0;
    sinfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sinfine_p_ = 0;
    cosfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    cosfine_p_ = 0;
    fail_ = false;
    message_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    message_ = '\0';
    trace_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    trace_ = '\0';
    // Initialize type pointers
    
  }

  virtual ~Fo_Dtwpgeometry_Master() = default;

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxpartials() const {
    return maxpartials_;
  }

  

  const int& maxfine() const {
    return maxfine_;
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

  

  const bool& do_upwelling() const {
    return do_upwelling_;
  }

  void do_upwelling(const bool& do_upwelling_in) {
    do_upwelling_ = do_upwelling_in;
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

  

  const bool& do_partials() const {
    return do_partials_;
  }

  void do_partials(const bool& do_partials_in) {
    do_partials_ = do_partials_in;
  }

  

  const int& ngeoms() const {
    return ngeoms_;
  }

  

  const int& nlayers() const {
    return nlayers_;
  }

  

  const int& npartials() const {
    return npartials_;
  }

  

  const int& nfine() const {
    return nfine_;
  }

  

  const blitz::Array<int, 1>& partial_layeridx() const {
    return partial_layeridx_;
  }

  void partial_layeridx(const blitz::Array<int, 1>& partial_layeridx_in) {
    partial_layeridx_ = partial_layeridx_in;
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

  

  const blitz::Array<double, 1>& partial_heights() const {
    return partial_heights_;
  }

  void partial_heights(const blitz::Array<double, 1>& partial_heights_in) {
    partial_heights_ = partial_heights_in;
  }

  

  const blitz::Array<double, 1>& mu1() const {
    return mu1_;
  }

  

  const blitz::Array<double, 1>& radii() const {
    return radii_;
  }

  void radii(const blitz::Array<double, 1>& radii_in) {
    radii_ = radii_in;
  }

  

  const blitz::Array<double, 2>& losw_paths() const {
    return losw_paths_;
  }

  

  const blitz::Array<double, 2>& alpha() const {
    return alpha_;
  }

  void alpha(const blitz::Array<double, 2>& alpha_in) {
    alpha_ = alpha_in;
  }

  

  const blitz::Array<double, 2>& sina() const {
    return sina_;
  }

  void sina(const blitz::Array<double, 2>& sina_in) {
    sina_ = sina_in;
  }

  

  const blitz::Array<double, 2>& cosa() const {
    return cosa_;
  }

  void cosa(const blitz::Array<double, 2>& cosa_in) {
    cosa_ = cosa_in;
  }

  

  const blitz::Array<double, 1>& radii_p() const {
    return radii_p_;
  }

  void radii_p(const blitz::Array<double, 1>& radii_p_in) {
    radii_p_ = radii_p_in;
  }

  

  const blitz::Array<double, 2>& losp_paths() const {
    return losp_paths_;
  }

  

  const blitz::Array<double, 2>& alpha_p() const {
    return alpha_p_;
  }

  void alpha_p(const blitz::Array<double, 2>& alpha_p_in) {
    alpha_p_ = alpha_p_in;
  }

  

  const blitz::Array<double, 2>& sina_p() const {
    return sina_p_;
  }

  void sina_p(const blitz::Array<double, 2>& sina_p_in) {
    sina_p_ = sina_p_in;
  }

  

  const blitz::Array<double, 2>& cosa_p() const {
    return cosa_p_;
  }

  void cosa_p(const blitz::Array<double, 2>& cosa_p_in) {
    cosa_p_ = cosa_p_in;
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

  

  const blitz::Array<double, 3>& radiifine() const {
    return radiifine_;
  }

  void radiifine(const blitz::Array<double, 3>& radiifine_in) {
    radiifine_ = radiifine_in;
  }

  

  const blitz::Array<double, 3>& alphafine() const {
    return alphafine_;
  }

  void alphafine(const blitz::Array<double, 3>& alphafine_in) {
    alphafine_ = alphafine_in;
  }

  

  const blitz::Array<double, 3>& sinfine() const {
    return sinfine_;
  }

  void sinfine(const blitz::Array<double, 3>& sinfine_in) {
    sinfine_ = sinfine_in;
  }

  

  const blitz::Array<double, 3>& cosfine() const {
    return cosfine_;
  }

  void cosfine(const blitz::Array<double, 3>& cosfine_in) {
    cosfine_ = cosfine_in;
  }

  

  const blitz::Array<int, 2>& nfinedivs_p() const {
    return nfinedivs_p_;
  }

  void nfinedivs_p(const blitz::Array<int, 2>& nfinedivs_p_in) {
    nfinedivs_p_ = nfinedivs_p_in;
  }

  

  const blitz::Array<double, 3>& xfine_p() const {
    return xfine_p_;
  }

  void xfine_p(const blitz::Array<double, 3>& xfine_p_in) {
    xfine_p_ = xfine_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_p() const {
    return wfine_p_;
  }

  void wfine_p(const blitz::Array<double, 3>& wfine_p_in) {
    wfine_p_ = wfine_p_in;
  }

  

  const blitz::Array<double, 3>& radiifine_p() const {
    return radiifine_p_;
  }

  void radiifine_p(const blitz::Array<double, 3>& radiifine_p_in) {
    radiifine_p_ = radiifine_p_in;
  }

  

  const blitz::Array<double, 3>& alphafine_p() const {
    return alphafine_p_;
  }

  void alphafine_p(const blitz::Array<double, 3>& alphafine_p_in) {
    alphafine_p_ = alphafine_p_in;
  }

  

  const blitz::Array<double, 3>& sinfine_p() const {
    return sinfine_p_;
  }

  void sinfine_p(const blitz::Array<double, 3>& sinfine_p_in) {
    sinfine_p_ = sinfine_p_in;
  }

  

  const blitz::Array<double, 3>& cosfine_p() const {
    return cosfine_p_;
  }

  void cosfine_p(const blitz::Array<double, 3>& cosfine_p_in) {
    cosfine_p_ = cosfine_p_in;
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
    
    fo_dtwpgeometry_master_m_fo_dtwpgeometry_master_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &dtr_, &eradius_, &do_upwelling_, &do_planpar_, &do_enhanced_ps_, &do_partials_, &ngeoms_, &nlayers_, &npartials_, &nfine_, partial_layeridx_.dataFirst(), heights_.dataFirst(), alpha_boa_.dataFirst(), partial_heights_.dataFirst(), mu1_.dataFirst(), radii_.dataFirst(), losw_paths_.dataFirst(), alpha_.dataFirst(), sina_.dataFirst(), cosa_.dataFirst(), radii_p_.dataFirst(), losp_paths_.dataFirst(), alpha_p_.dataFirst(), sina_p_.dataFirst(), cosa_p_.dataFirst(), nfinedivs_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), radiifine_.dataFirst(), alphafine_.dataFirst(), sinfine_.dataFirst(), cosfine_.dataFirst(), nfinedivs_p_.dataFirst(), xfine_p_.dataFirst(), wfine_p_.dataFirst(), radiifine_p_.dataFirst(), alphafine_p_.dataFirst(), sinfine_p_.dataFirst(), cosfine_p_.dataFirst(), &fail_, &message_len, message_.dataFirst(), &trace_len, trace_.dataFirst());
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "Fo_Dtwpgeometry_Master:" << std::endl
      << "        maxgeoms: " << maxgeoms()  << std::endl
      << "       maxlayers: " << maxlayers()  << std::endl
      << "     maxpartials: " << maxpartials()  << std::endl
      << "         maxfine: " << maxfine()  << std::endl
      << "             dtr: " << dtr()  << std::endl
      << "         eradius: " << eradius()  << std::endl
      << "    do_upwelling: " << do_upwelling()  << std::endl
      << "      do_planpar: " << do_planpar()  << std::endl
      << "  do_enhanced_ps: " << do_enhanced_ps()  << std::endl
      << "     do_partials: " << do_partials()  << std::endl
      << "          ngeoms: " << ngeoms()  << std::endl
      << "         nlayers: " << nlayers()  << std::endl
      << "       npartials: " << npartials()  << std::endl
      << "           nfine: " << nfine()  << std::endl
      << "partial_layeridx: " << std::endl << partial_layeridx()  << std::endl
      << "         heights: " << std::endl << heights()  << std::endl
      << "       alpha_boa: " << std::endl << alpha_boa()  << std::endl
      << " partial_heights: " << std::endl << partial_heights()  << std::endl
      << "             mu1: " << std::endl << mu1()  << std::endl
      << "           radii: " << std::endl << radii()  << std::endl
      << "      losw_paths: " << std::endl << losw_paths()  << std::endl
      << "           alpha: " << std::endl << alpha()  << std::endl
      << "            sina: " << std::endl << sina()  << std::endl
      << "            cosa: " << std::endl << cosa()  << std::endl
      << "         radii_p: " << std::endl << radii_p()  << std::endl
      << "      losp_paths: " << std::endl << losp_paths()  << std::endl
      << "         alpha_p: " << std::endl << alpha_p()  << std::endl
      << "          sina_p: " << std::endl << sina_p()  << std::endl
      << "          cosa_p: " << std::endl << cosa_p()  << std::endl
      << "       nfinedivs: " << std::endl << nfinedivs()  << std::endl
      << "           xfine: " << std::endl << xfine()  << std::endl
      << "           wfine: " << std::endl << wfine()  << std::endl
      << "       radiifine: " << std::endl << radiifine()  << std::endl
      << "       alphafine: " << std::endl << alphafine()  << std::endl
      << "         sinfine: " << std::endl << sinfine()  << std::endl
      << "         cosfine: " << std::endl << cosfine()  << std::endl
      << "     nfinedivs_p: " << std::endl << nfinedivs_p()  << std::endl
      << "         xfine_p: " << std::endl << xfine_p()  << std::endl
      << "         wfine_p: " << std::endl << wfine_p()  << std::endl
      << "     radiifine_p: " << std::endl << radiifine_p()  << std::endl
      << "     alphafine_p: " << std::endl << alphafine_p()  << std::endl
      << "       sinfine_p: " << std::endl << sinfine_p()  << std::endl
      << "       cosfine_p: " << std::endl << cosfine_p()  << std::endl
      << "            fail: " << fail()  << std::endl
      << "         message: " << "\"" << message() << "\"" << std::endl
      << "           trace: " << "\"" << trace() << "\"" << std::endl;

  }

private:
  int maxgeoms_;
  int maxlayers_;
  int maxpartials_;
  int maxfine_;
  double dtr_;
  double eradius_;
  bool do_upwelling_;
  bool do_planpar_;
  bool do_enhanced_ps_;
  bool do_partials_;
  int ngeoms_;
  int nlayers_;
  int npartials_;
  int nfine_;
  blitz::Array<int, 1> partial_layeridx_;
  blitz::Array<double, 1> heights_;
  blitz::Array<double, 1> alpha_boa_;
  blitz::Array<double, 1> partial_heights_;
  blitz::Array<double, 1> mu1_;
  blitz::Array<double, 1> radii_;
  blitz::Array<double, 2> losw_paths_;
  blitz::Array<double, 2> alpha_;
  blitz::Array<double, 2> sina_;
  blitz::Array<double, 2> cosa_;
  blitz::Array<double, 1> radii_p_;
  blitz::Array<double, 2> losp_paths_;
  blitz::Array<double, 2> alpha_p_;
  blitz::Array<double, 2> sina_p_;
  blitz::Array<double, 2> cosa_p_;
  blitz::Array<int, 2> nfinedivs_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> radiifine_;
  blitz::Array<double, 3> alphafine_;
  blitz::Array<double, 3> sinfine_;
  blitz::Array<double, 3> cosfine_;
  blitz::Array<int, 2> nfinedivs_p_;
  blitz::Array<double, 3> xfine_p_;
  blitz::Array<double, 3> wfine_p_;
  blitz::Array<double, 3> radiifine_p_;
  blitz::Array<double, 3> alphafine_p_;
  blitz::Array<double, 3> sinfine_p_;
  blitz::Array<double, 3> cosfine_p_;
  bool fail_;
  blitz::Array<char, 1> message_;
  blitz::Array<char, 1> trace_;

  // Serialization support
  Fo_Dtwpgeometry_Master() {}
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "fo_sswpgeometry_master_m" in file: "FO_SSWPgeometry_master.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_sswpgeometry_master_m_fo_sswpgeometry_master_wrap(const int* maxgeoms_in, const int* maxszas_in, const int* maxvzas_in, const int* maxazms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const bool* do_obsgeom_in, const bool* do_doublet_in, const bool* do_chapman_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const bool* do_partials_in, const int* ngeoms_in, const int* nszas_in, const int* nvzas_in, const int* nazms_in, const int* nlayers_in, const int* nfine_in, const int* npartials_in, const int* partial_layeridx_in, const double* dtr_in, const double* pie_in, const double* vsign_in, const double* eradius_in, const int* nv_offset_in, const int* na_offset_in, const int* nd_offset_in, const double* heights_in, const double* partial_heights_in, const double* obsgeom_boa_in, const double* alpha_boa_in, const double* theta_boa_in, const double* phi_boa_in, const bool* donadir_in, const double* raycon_in, const double* mu0_in, const double* mu1_in, const double* cosscat_in, const double* radii_in, const double* losw_paths_in, const double* alpha_in, const double* sina_in, const double* cosa_in, const double* sunpaths_in, const int* ntraverse_in, const double* chapfacs_in, const double* theta_all_in, const double* radii_p_in, const double* losp_paths_in, const double* alpha_p_in, const double* sina_p_in, const double* cosa_p_in, const double* sunpaths_p_in, const int* ntraverse_p_in, const double* chapfacs_p_in, const int* nfinedivs_in, const double* xfine_in, const double* wfine_in, const double* radiifine_in, const double* alphafine_in, const double* sinfine_in, const double* cosfine_in, const double* sunpathsfine_in, const int* ntraversefine_in, const int* nfinedivs_p_in, const double* xfine_p_in, const double* wfine_p_in, const double* radiifine_p_in, const double* alphafine_p_in, const double* sinfine_p_in, const double* cosfine_p_in, const double* sunpathsfine_p_in, const int* ntraversefine_p_in, const bool* fail_in, const int* message_len, const char* message_in, const int* trace_len, const char* trace_in);
}

class Fo_Sswpgeometry_Master : public Printable<Fo_Sswpgeometry_Master> {

public:
  Fo_Sswpgeometry_Master(const int& maxgeoms_in, const int& maxszas_in, const int& maxvzas_in, const int& maxazms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& ngeoms_in, const int& nszas_in, const int& nvzas_in, const int& nazms_in, const int& nlayers_in, const int& nfine_in, const int& npartials_in) : maxgeoms_(maxgeoms_in), maxszas_(maxszas_in), maxvzas_(maxvzas_in), maxazms_(maxazms_in), maxlayers_(maxlayers_in), maxpartials_(maxpartials_in), maxfine_(maxfine_in), ngeoms_(ngeoms_in), nszas_(nszas_in), nvzas_(nvzas_in), nazms_(nazms_in), nlayers_(nlayers_in), nfine_(nfine_in), npartials_(npartials_in) 
  { 
    do_obsgeom_ = false;
    do_doublet_ = false;
    do_chapman_ = false;
    do_planpar_ = false;
    do_enhanced_ps_ = false;
    do_partials_ = false;
    partial_layeridx_.reference( blitz::Array<int, 1>(maxpartials_, blitz::ColumnMajorArray<1>()) );
    partial_layeridx_ = 0;
    dtr_ = 0;
    pie_ = 0;
    vsign_ = 0;
    eradius_ = 0;
    nv_offset_.reference( blitz::Array<int, 1>(maxszas_, blitz::ColumnMajorArray<1>()) );
    nv_offset_ = 0;
    na_offset_.reference( blitz::Array<int, 2>(maxszas_, maxvzas_, blitz::ColumnMajorArray<2>()) );
    na_offset_ = 0;
    nd_offset_.reference( blitz::Array<int, 1>(maxszas_, blitz::ColumnMajorArray<1>()) );
    nd_offset_ = 0;
    heights_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    heights_ = 0;
    partial_heights_.reference( blitz::Array<double, 1>(maxpartials_, blitz::ColumnMajorArray<1>()) );
    partial_heights_ = 0;
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
    raycon_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    raycon_ = 0;
    mu0_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu0_ = 0;
    mu1_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu1_ = 0;
    cosscat_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    cosscat_ = 0;
    radii_.reference( blitz::Array<double, 1>(maxlayers_-0+1, blitz::ColumnMajorArray<1>()) );
    radii_ = 0;
    losw_paths_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losw_paths_ = 0;
    alpha_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    alpha_ = 0;
    sina_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    sina_ = 0;
    cosa_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cosa_ = 0;
    sunpaths_.reference( blitz::Array<double, 3>(maxlayers_-0+1, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_ = 0;
    ntraverse_.reference( blitz::Array<int, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_ = 0;
    chapfacs_.reference( blitz::Array<double, 3>(maxlayers_, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    chapfacs_ = 0;
    theta_all_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    theta_all_ = 0;
    radii_p_.reference( blitz::Array<double, 1>(maxpartials_, blitz::ColumnMajorArray<1>()) );
    radii_p_ = 0;
    losp_paths_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losp_paths_ = 0;
    alpha_p_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    alpha_p_ = 0;
    sina_p_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    sina_p_ = 0;
    cosa_p_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cosa_p_ = 0;
    sunpaths_p_.reference( blitz::Array<double, 3>(maxpartials_, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_p_ = 0;
    ntraverse_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_p_ = 0;
    chapfacs_p_.reference( blitz::Array<double, 3>(maxpartials_, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    chapfacs_p_ = 0;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    radiifine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    radiifine_ = 0;
    alphafine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    alphafine_ = 0;
    sinfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sinfine_ = 0;
    cosfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    cosfine_ = 0;
    sunpathsfine_.reference( blitz::Array<double, 4>(maxlayers_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_ = 0;
    ntraversefine_.reference( blitz::Array<int, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_ = 0;
    nfinedivs_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_p_ = 0;
    xfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_p_ = 0;
    wfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_p_ = 0;
    radiifine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    radiifine_p_ = 0;
    alphafine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    alphafine_p_ = 0;
    sinfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sinfine_p_ = 0;
    cosfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    cosfine_p_ = 0;
    sunpathsfine_p_.reference( blitz::Array<double, 4>(maxpartials_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_p_ = 0;
    ntraversefine_p_.reference( blitz::Array<int, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_p_ = 0;
    fail_ = false;
    message_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    message_ = '\0';
    trace_.reference( blitz::Array<char, 1>(101, blitz::ColumnMajorArray<1>()) );
    trace_ = '\0';
    // Initialize type pointers
    
  }

  virtual ~Fo_Sswpgeometry_Master() = default;

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

  

  const int& maxpartials() const {
    return maxpartials_;
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

  

  const bool& do_doublet() const {
    return do_doublet_;
  }

  void do_doublet(const bool& do_doublet_in) {
    do_doublet_ = do_doublet_in;
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

  

  const bool& do_partials() const {
    return do_partials_;
  }

  void do_partials(const bool& do_partials_in) {
    do_partials_ = do_partials_in;
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

  

  const int& npartials() const {
    return npartials_;
  }

  

  const blitz::Array<int, 1>& partial_layeridx() const {
    return partial_layeridx_;
  }

  void partial_layeridx(const blitz::Array<int, 1>& partial_layeridx_in) {
    partial_layeridx_ = partial_layeridx_in;
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

  

  const blitz::Array<int, 1>& nv_offset() const {
    return nv_offset_;
  }

  void nv_offset(const blitz::Array<int, 1>& nv_offset_in) {
    nv_offset_ = nv_offset_in;
  }

  

  const blitz::Array<int, 2>& na_offset() const {
    return na_offset_;
  }

  void na_offset(const blitz::Array<int, 2>& na_offset_in) {
    na_offset_ = na_offset_in;
  }

  

  const blitz::Array<int, 1>& nd_offset() const {
    return nd_offset_;
  }

  void nd_offset(const blitz::Array<int, 1>& nd_offset_in) {
    nd_offset_ = nd_offset_in;
  }

  

  const blitz::Array<double, 1>& heights() const {
    return heights_;
  }

  void heights(const blitz::Array<double, 1>& heights_in) {
    heights_ = heights_in;
  }

  

  const blitz::Array<double, 1>& partial_heights() const {
    return partial_heights_;
  }

  void partial_heights(const blitz::Array<double, 1>& partial_heights_in) {
    partial_heights_ = partial_heights_in;
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

  

  const blitz::Array<double, 1>& raycon() const {
    return raycon_;
  }

  void raycon(const blitz::Array<double, 1>& raycon_in) {
    raycon_ = raycon_in;
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

  

  const blitz::Array<double, 1>& radii() const {
    return radii_;
  }

  void radii(const blitz::Array<double, 1>& radii_in) {
    radii_ = radii_in;
  }

  

  const blitz::Array<double, 2>& losw_paths() const {
    return losw_paths_;
  }

  

  const blitz::Array<double, 2>& alpha() const {
    return alpha_;
  }

  void alpha(const blitz::Array<double, 2>& alpha_in) {
    alpha_ = alpha_in;
  }

  

  const blitz::Array<double, 2>& sina() const {
    return sina_;
  }

  void sina(const blitz::Array<double, 2>& sina_in) {
    sina_ = sina_in;
  }

  

  const blitz::Array<double, 2>& cosa() const {
    return cosa_;
  }

  void cosa(const blitz::Array<double, 2>& cosa_in) {
    cosa_ = cosa_in;
  }

  

  const blitz::Array<double, 3>& sunpaths() const {
    return sunpaths_;
  }

  

  const blitz::Array<int, 2>& ntraverse() const {
    return ntraverse_;
  }

  

  const blitz::Array<double, 3>& chapfacs() const {
    return chapfacs_;
  }

  

  const blitz::Array<double, 2>& theta_all() const {
    return theta_all_;
  }

  

  const blitz::Array<double, 1>& radii_p() const {
    return radii_p_;
  }

  void radii_p(const blitz::Array<double, 1>& radii_p_in) {
    radii_p_ = radii_p_in;
  }

  

  const blitz::Array<double, 2>& losp_paths() const {
    return losp_paths_;
  }

  

  const blitz::Array<double, 2>& alpha_p() const {
    return alpha_p_;
  }

  void alpha_p(const blitz::Array<double, 2>& alpha_p_in) {
    alpha_p_ = alpha_p_in;
  }

  

  const blitz::Array<double, 2>& sina_p() const {
    return sina_p_;
  }

  void sina_p(const blitz::Array<double, 2>& sina_p_in) {
    sina_p_ = sina_p_in;
  }

  

  const blitz::Array<double, 2>& cosa_p() const {
    return cosa_p_;
  }

  void cosa_p(const blitz::Array<double, 2>& cosa_p_in) {
    cosa_p_ = cosa_p_in;
  }

  

  const blitz::Array<double, 3>& sunpaths_p() const {
    return sunpaths_p_;
  }

  

  const blitz::Array<int, 2>& ntraverse_p() const {
    return ntraverse_p_;
  }

  

  const blitz::Array<double, 3>& chapfacs_p() const {
    return chapfacs_p_;
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

  

  const blitz::Array<double, 3>& radiifine() const {
    return radiifine_;
  }

  void radiifine(const blitz::Array<double, 3>& radiifine_in) {
    radiifine_ = radiifine_in;
  }

  

  const blitz::Array<double, 3>& alphafine() const {
    return alphafine_;
  }

  void alphafine(const blitz::Array<double, 3>& alphafine_in) {
    alphafine_ = alphafine_in;
  }

  

  const blitz::Array<double, 3>& sinfine() const {
    return sinfine_;
  }

  void sinfine(const blitz::Array<double, 3>& sinfine_in) {
    sinfine_ = sinfine_in;
  }

  

  const blitz::Array<double, 3>& cosfine() const {
    return cosfine_;
  }

  void cosfine(const blitz::Array<double, 3>& cosfine_in) {
    cosfine_ = cosfine_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine() const {
    return sunpathsfine_;
  }

  

  const blitz::Array<int, 3>& ntraversefine() const {
    return ntraversefine_;
  }

  

  const blitz::Array<int, 2>& nfinedivs_p() const {
    return nfinedivs_p_;
  }

  void nfinedivs_p(const blitz::Array<int, 2>& nfinedivs_p_in) {
    nfinedivs_p_ = nfinedivs_p_in;
  }

  

  const blitz::Array<double, 3>& xfine_p() const {
    return xfine_p_;
  }

  void xfine_p(const blitz::Array<double, 3>& xfine_p_in) {
    xfine_p_ = xfine_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_p() const {
    return wfine_p_;
  }

  void wfine_p(const blitz::Array<double, 3>& wfine_p_in) {
    wfine_p_ = wfine_p_in;
  }

  

  const blitz::Array<double, 3>& radiifine_p() const {
    return radiifine_p_;
  }

  void radiifine_p(const blitz::Array<double, 3>& radiifine_p_in) {
    radiifine_p_ = radiifine_p_in;
  }

  

  const blitz::Array<double, 3>& alphafine_p() const {
    return alphafine_p_;
  }

  void alphafine_p(const blitz::Array<double, 3>& alphafine_p_in) {
    alphafine_p_ = alphafine_p_in;
  }

  

  const blitz::Array<double, 3>& sinfine_p() const {
    return sinfine_p_;
  }

  void sinfine_p(const blitz::Array<double, 3>& sinfine_p_in) {
    sinfine_p_ = sinfine_p_in;
  }

  

  const blitz::Array<double, 3>& cosfine_p() const {
    return cosfine_p_;
  }

  void cosfine_p(const blitz::Array<double, 3>& cosfine_p_in) {
    cosfine_p_ = cosfine_p_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine_p() const {
    return sunpathsfine_p_;
  }

  

  const blitz::Array<int, 3>& ntraversefine_p() const {
    return ntraversefine_p_;
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
    
    fo_sswpgeometry_master_m_fo_sswpgeometry_master_wrap(&maxgeoms_, &maxszas_, &maxvzas_, &maxazms_, &maxlayers_, &maxpartials_, &maxfine_, &do_obsgeom_, &do_doublet_, &do_chapman_, &do_planpar_, &do_enhanced_ps_, &do_partials_, &ngeoms_, &nszas_, &nvzas_, &nazms_, &nlayers_, &nfine_, &npartials_, partial_layeridx_.dataFirst(), &dtr_, &pie_, &vsign_, &eradius_, nv_offset_.dataFirst(), na_offset_.dataFirst(), nd_offset_.dataFirst(), heights_.dataFirst(), partial_heights_.dataFirst(), obsgeom_boa_.dataFirst(), alpha_boa_.dataFirst(), theta_boa_.dataFirst(), phi_boa_.dataFirst(), donadir_.dataFirst(), raycon_.dataFirst(), mu0_.dataFirst(), mu1_.dataFirst(), cosscat_.dataFirst(), radii_.dataFirst(), losw_paths_.dataFirst(), alpha_.dataFirst(), sina_.dataFirst(), cosa_.dataFirst(), sunpaths_.dataFirst(), ntraverse_.dataFirst(), chapfacs_.dataFirst(), theta_all_.dataFirst(), radii_p_.dataFirst(), losp_paths_.dataFirst(), alpha_p_.dataFirst(), sina_p_.dataFirst(), cosa_p_.dataFirst(), sunpaths_p_.dataFirst(), ntraverse_p_.dataFirst(), chapfacs_p_.dataFirst(), nfinedivs_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), radiifine_.dataFirst(), alphafine_.dataFirst(), sinfine_.dataFirst(), cosfine_.dataFirst(), sunpathsfine_.dataFirst(), ntraversefine_.dataFirst(), nfinedivs_p_.dataFirst(), xfine_p_.dataFirst(), wfine_p_.dataFirst(), radiifine_p_.dataFirst(), alphafine_p_.dataFirst(), sinfine_p_.dataFirst(), cosfine_p_.dataFirst(), sunpathsfine_p_.dataFirst(), ntraversefine_p_.dataFirst(), &fail_, &message_len, message_.dataFirst(), &trace_len, trace_.dataFirst());
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "Fo_Sswpgeometry_Master:" << std::endl
      << "        maxgeoms: " << maxgeoms()  << std::endl
      << "         maxszas: " << maxszas()  << std::endl
      << "         maxvzas: " << maxvzas()  << std::endl
      << "         maxazms: " << maxazms()  << std::endl
      << "       maxlayers: " << maxlayers()  << std::endl
      << "     maxpartials: " << maxpartials()  << std::endl
      << "         maxfine: " << maxfine()  << std::endl
      << "      do_obsgeom: " << do_obsgeom()  << std::endl
      << "      do_doublet: " << do_doublet()  << std::endl
      << "      do_chapman: " << do_chapman()  << std::endl
      << "      do_planpar: " << do_planpar()  << std::endl
      << "  do_enhanced_ps: " << do_enhanced_ps()  << std::endl
      << "     do_partials: " << do_partials()  << std::endl
      << "          ngeoms: " << ngeoms()  << std::endl
      << "           nszas: " << nszas()  << std::endl
      << "           nvzas: " << nvzas()  << std::endl
      << "           nazms: " << nazms()  << std::endl
      << "         nlayers: " << nlayers()  << std::endl
      << "           nfine: " << nfine()  << std::endl
      << "       npartials: " << npartials()  << std::endl
      << "partial_layeridx: " << std::endl << partial_layeridx()  << std::endl
      << "             dtr: " << dtr()  << std::endl
      << "             pie: " << pie()  << std::endl
      << "           vsign: " << vsign()  << std::endl
      << "         eradius: " << eradius()  << std::endl
      << "       nv_offset: " << std::endl << nv_offset()  << std::endl
      << "       na_offset: " << std::endl << na_offset()  << std::endl
      << "       nd_offset: " << std::endl << nd_offset()  << std::endl
      << "         heights: " << std::endl << heights()  << std::endl
      << " partial_heights: " << std::endl << partial_heights()  << std::endl
      << "     obsgeom_boa: " << std::endl << obsgeom_boa()  << std::endl
      << "       alpha_boa: " << std::endl << alpha_boa()  << std::endl
      << "       theta_boa: " << std::endl << theta_boa()  << std::endl
      << "         phi_boa: " << std::endl << phi_boa()  << std::endl
      << "         donadir: " << std::endl << donadir()  << std::endl
      << "          raycon: " << std::endl << raycon()  << std::endl
      << "             mu0: " << std::endl << mu0()  << std::endl
      << "             mu1: " << std::endl << mu1()  << std::endl
      << "         cosscat: " << std::endl << cosscat()  << std::endl
      << "           radii: " << std::endl << radii()  << std::endl
      << "      losw_paths: " << std::endl << losw_paths()  << std::endl
      << "           alpha: " << std::endl << alpha()  << std::endl
      << "            sina: " << std::endl << sina()  << std::endl
      << "            cosa: " << std::endl << cosa()  << std::endl
      << "        sunpaths: " << std::endl << sunpaths()  << std::endl
      << "       ntraverse: " << std::endl << ntraverse()  << std::endl
      << "        chapfacs: " << std::endl << chapfacs()  << std::endl
      << "       theta_all: " << std::endl << theta_all()  << std::endl
      << "         radii_p: " << std::endl << radii_p()  << std::endl
      << "      losp_paths: " << std::endl << losp_paths()  << std::endl
      << "         alpha_p: " << std::endl << alpha_p()  << std::endl
      << "          sina_p: " << std::endl << sina_p()  << std::endl
      << "          cosa_p: " << std::endl << cosa_p()  << std::endl
      << "      sunpaths_p: " << std::endl << sunpaths_p()  << std::endl
      << "     ntraverse_p: " << std::endl << ntraverse_p()  << std::endl
      << "      chapfacs_p: " << std::endl << chapfacs_p()  << std::endl
      << "       nfinedivs: " << std::endl << nfinedivs()  << std::endl
      << "           xfine: " << std::endl << xfine()  << std::endl
      << "           wfine: " << std::endl << wfine()  << std::endl
      << "       radiifine: " << std::endl << radiifine()  << std::endl
      << "       alphafine: " << std::endl << alphafine()  << std::endl
      << "         sinfine: " << std::endl << sinfine()  << std::endl
      << "         cosfine: " << std::endl << cosfine()  << std::endl
      << "    sunpathsfine: " << std::endl << sunpathsfine()  << std::endl
      << "   ntraversefine: " << std::endl << ntraversefine()  << std::endl
      << "     nfinedivs_p: " << std::endl << nfinedivs_p()  << std::endl
      << "         xfine_p: " << std::endl << xfine_p()  << std::endl
      << "         wfine_p: " << std::endl << wfine_p()  << std::endl
      << "     radiifine_p: " << std::endl << radiifine_p()  << std::endl
      << "     alphafine_p: " << std::endl << alphafine_p()  << std::endl
      << "       sinfine_p: " << std::endl << sinfine_p()  << std::endl
      << "       cosfine_p: " << std::endl << cosfine_p()  << std::endl
      << "  sunpathsfine_p: " << std::endl << sunpathsfine_p()  << std::endl
      << " ntraversefine_p: " << std::endl << ntraversefine_p()  << std::endl
      << "            fail: " << fail()  << std::endl
      << "         message: " << "\"" << message() << "\"" << std::endl
      << "           trace: " << "\"" << trace() << "\"" << std::endl;

  }

private:
  int maxgeoms_;
  int maxszas_;
  int maxvzas_;
  int maxazms_;
  int maxlayers_;
  int maxpartials_;
  int maxfine_;
  bool do_obsgeom_;
  bool do_doublet_;
  bool do_chapman_;
  bool do_planpar_;
  bool do_enhanced_ps_;
  bool do_partials_;
  int ngeoms_;
  int nszas_;
  int nvzas_;
  int nazms_;
  int nlayers_;
  int nfine_;
  int npartials_;
  blitz::Array<int, 1> partial_layeridx_;
  double dtr_;
  double pie_;
  double vsign_;
  double eradius_;
  blitz::Array<int, 1> nv_offset_;
  blitz::Array<int, 2> na_offset_;
  blitz::Array<int, 1> nd_offset_;
  blitz::Array<double, 1> heights_;
  blitz::Array<double, 1> partial_heights_;
  blitz::Array<double, 2> obsgeom_boa_;
  blitz::Array<double, 1> alpha_boa_;
  blitz::Array<double, 1> theta_boa_;
  blitz::Array<double, 1> phi_boa_;
  blitz::Array<bool, 1> donadir_;
  blitz::Array<double, 1> raycon_;
  blitz::Array<double, 1> mu0_;
  blitz::Array<double, 1> mu1_;
  blitz::Array<double, 1> cosscat_;
  blitz::Array<double, 1> radii_;
  blitz::Array<double, 2> losw_paths_;
  blitz::Array<double, 2> alpha_;
  blitz::Array<double, 2> sina_;
  blitz::Array<double, 2> cosa_;
  blitz::Array<double, 3> sunpaths_;
  blitz::Array<int, 2> ntraverse_;
  blitz::Array<double, 3> chapfacs_;
  blitz::Array<double, 2> theta_all_;
  blitz::Array<double, 1> radii_p_;
  blitz::Array<double, 2> losp_paths_;
  blitz::Array<double, 2> alpha_p_;
  blitz::Array<double, 2> sina_p_;
  blitz::Array<double, 2> cosa_p_;
  blitz::Array<double, 3> sunpaths_p_;
  blitz::Array<int, 2> ntraverse_p_;
  blitz::Array<double, 3> chapfacs_p_;
  blitz::Array<int, 2> nfinedivs_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> radiifine_;
  blitz::Array<double, 3> alphafine_;
  blitz::Array<double, 3> sinfine_;
  blitz::Array<double, 3> cosfine_;
  blitz::Array<double, 4> sunpathsfine_;
  blitz::Array<int, 3> ntraversefine_;
  blitz::Array<int, 2> nfinedivs_p_;
  blitz::Array<double, 3> xfine_p_;
  blitz::Array<double, 3> wfine_p_;
  blitz::Array<double, 3> radiifine_p_;
  blitz::Array<double, 3> alphafine_p_;
  blitz::Array<double, 3> sinfine_p_;
  blitz::Array<double, 3> cosfine_p_;
  blitz::Array<double, 4> sunpathsfine_p_;
  blitz::Array<int, 3> ntraversefine_p_;
  bool fail_;
  blitz::Array<char, 1> message_;
  blitz::Array<char, 1> trace_;

  // Serialization support
  Fo_Sswpgeometry_Master() {}
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "fo_scalarss_rtcalcs_i_m" in file: "FO_ScalarSS_RTCalcs_I.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_scalarss_rtcalcs_i_m_ss_integral_i_dn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* maxmoments_input_in, const int* max_user_levels_in, const bool* do_deltam_scaling_in, const bool* do_phasfunc_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const double* flux_in, const bool* do_sources_dn_in, const bool* do_sources_dn_p_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* nmoments_input_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* phasmoms_in, const double* phasfunc_dn_in, const double* mu1_in, const double* legpoly_dn_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_in, const double* wfine_in, const double* sunpaths_in, const int* ntraverse_in, const double* sunpathsfine_in, const int* ntraversefine_in, const double* xfine_p_in, const double* wfine_p_in, const double* sunpaths_p_in, const int* ntraverse_p_in, const double* sunpathsfine_p_in, const int* ntraversefine_p_in, const double* intensity_dn_in, const double* cumsource_dn_in, const double* lostrans_dn_in);
  void fo_scalarss_rtcalcs_i_m_ss_integral_i_up_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* maxmoments_input_in, const int* max_user_levels_in, const bool* do_deltam_scaling_in, const bool* do_phasfunc_in, const bool* do_surface_leaving_in, const bool* do_water_leaving_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const double* flux_in, const bool* do_sources_up_in, const bool* do_sources_up_p_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* nmoments_input_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* phasmoms_in, const double* phasfunc_up_in, const double* reflec_in, const double* slterm_in, const double* mu0_in, const double* mu1_in, const double* legpoly_up_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_in, const double* wfine_in, const double* sunpaths_in, const int* ntraverse_in, const double* sunpathsfine_in, const int* ntraversefine_in, const double* xfine_p_in, const double* wfine_p_in, const double* sunpaths_p_in, const int* ntraverse_p_in, const double* sunpathsfine_p_in, const int* ntraversefine_p_in, const double* intensity_up_in, const double* intensity_db_in, const double* cumsource_up_in, const double* cumtrans_in, const double* lostrans_up_in);
  void fo_scalarss_rtcalcs_i_m_ss_integral_i_updn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* maxmoments_input_in, const int* max_user_levels_in, const bool* do_upwelling_in, const bool* do_dnwelling_in, const bool* do_deltam_scaling_in, const bool* do_phasfunc_in, const bool* do_surface_leaving_in, const bool* do_water_leaving_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const bool* do_sources_up_in, const bool* do_sources_up_p_in, const bool* do_sources_dn_in, const bool* do_sources_dn_p_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* nmoments_input_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* flux_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* phasmoms_in, const double* phasfunc_up_in, const double* phasfunc_dn_in, const double* reflec_in, const double* slterm_in, const double* mu0_in, const double* mu1_in, const double* legpoly_up_in, const double* legpoly_dn_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_up_in, const double* wfine_up_in, const double* sunpaths_up_in, const int* ntraverse_up_in, const double* sunpathsfine_up_in, const int* ntraversefine_up_in, const double* xfine_dn_in, const double* wfine_dn_in, const double* sunpaths_dn_in, const int* ntraverse_dn_in, const double* sunpathsfine_dn_in, const int* ntraversefine_dn_in, const double* xfine_up_p_in, const double* wfine_up_p_in, const double* sunpaths_up_p_in, const int* ntraverse_up_p_in, const double* sunpathsfine_up_p_in, const int* ntraversefine_up_p_in, const double* xfine_dn_p_in, const double* wfine_dn_p_in, const double* sunpaths_dn_p_in, const int* ntraverse_dn_p_in, const double* sunpathsfine_dn_p_in, const int* ntraversefine_dn_p_in, const double* intensity_up_in, const double* intensity_db_in, const double* cumsource_up_in, const double* cumtrans_in, const double* lostrans_up_in, const double* intensity_dn_in, const double* cumsource_dn_in, const double* lostrans_dn_in);
}

class Fo_Scalarss_Rtcalcs_I : public Printable<Fo_Scalarss_Rtcalcs_I> {

public:
  Fo_Scalarss_Rtcalcs_I(const int& maxgeoms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& maxmoments_input_in, const int& max_user_levels_in, const int& ngeoms_in, const int& nlayers_in, const int& nmoments_input_in, const int& n_user_levels_in, const int& npartials_in) : maxgeoms_(maxgeoms_in), maxlayers_(maxlayers_in), maxpartials_(maxpartials_in), maxfine_(maxfine_in), maxmoments_input_(maxmoments_input_in), max_user_levels_(max_user_levels_in), ngeoms_(ngeoms_in), nlayers_(nlayers_in), nmoments_input_(nmoments_input_in), n_user_levels_(n_user_levels_in), npartials_(npartials_in) 
  { 
    do_deltam_scaling_ = false;
    do_phasfunc_ = false;
    do_partials_ = false;
    do_planpar_ = false;
    do_enhanced_ps_ = false;
    flux_ = 0;
    do_sources_dn_.reference( blitz::Array<bool, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_dn_ = false;
    do_sources_dn_p_.reference( blitz::Array<bool, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_dn_p_ = false;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    user_levels_.reference( blitz::Array<int, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    user_levels_ = 0;
    nfinedivs_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_p_ = 0;
    partial_outindex_.reference( blitz::Array<int, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    partial_outindex_ = 0;
    partial_outflag_.reference( blitz::Array<bool, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    partial_outflag_ = false;
    partial_layeridx_.reference( blitz::Array<int, 1>(maxpartials_, blitz::ColumnMajorArray<1>()) );
    partial_layeridx_ = 0;
    extinction_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    extinction_ = 0;
    deltaus_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    deltaus_ = 0;
    omega_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    omega_ = 0;
    truncfac_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    truncfac_ = 0;
    phasmoms_.reference( blitz::Array<double, 2>(maxlayers_, maxmoments_input_-0+1, blitz::ColumnMajorArray<2>()) );
    phasmoms_ = 0;
    phasfunc_dn_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    phasfunc_dn_ = 0;
    mu1_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu1_ = 0;
    legpoly_dn_.reference( blitz::Array<double, 2>(maxmoments_input_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    legpoly_dn_ = 0;
    losw_paths_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losw_paths_ = 0;
    losp_paths_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losp_paths_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    sunpaths_.reference( blitz::Array<double, 3>(maxlayers_-0+1, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_ = 0;
    ntraverse_.reference( blitz::Array<int, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_ = 0;
    sunpathsfine_.reference( blitz::Array<double, 4>(maxlayers_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_ = 0;
    ntraversefine_.reference( blitz::Array<int, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_ = 0;
    xfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_p_ = 0;
    wfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_p_ = 0;
    sunpaths_p_.reference( blitz::Array<double, 3>(maxpartials_, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_p_ = 0;
    ntraverse_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_p_ = 0;
    sunpathsfine_p_.reference( blitz::Array<double, 4>(maxpartials_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_p_ = 0;
    ntraversefine_p_.reference( blitz::Array<int, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_p_ = 0;
    intensity_dn_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dn_ = 0;
    cumsource_dn_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cumsource_dn_ = 0;
    lostrans_dn_.reference( blitz::Array<double, 2>(maxgeoms_, maxlayers_, blitz::ColumnMajorArray<2>()) );
    lostrans_dn_ = 0;
    do_surface_leaving_ = false;
    do_water_leaving_ = false;
    do_sources_up_.reference( blitz::Array<bool, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_up_ = false;
    do_sources_up_p_.reference( blitz::Array<bool, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_up_p_ = false;
    phasfunc_up_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    phasfunc_up_ = 0;
    reflec_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    reflec_ = 0;
    slterm_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    slterm_ = 0;
    mu0_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu0_ = 0;
    legpoly_up_.reference( blitz::Array<double, 2>(maxmoments_input_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    legpoly_up_ = 0;
    intensity_up_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_up_ = 0;
    intensity_db_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_db_ = 0;
    cumsource_up_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cumsource_up_ = 0;
    cumtrans_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cumtrans_ = 0;
    lostrans_up_.reference( blitz::Array<double, 2>(maxgeoms_, maxlayers_, blitz::ColumnMajorArray<2>()) );
    lostrans_up_ = 0;
    do_upwelling_ = false;
    do_dnwelling_ = false;
    xfine_up_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_up_ = 0;
    wfine_up_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_up_ = 0;
    sunpaths_up_.reference( blitz::Array<double, 3>(maxlayers_-0+1, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_up_ = 0;
    ntraverse_up_.reference( blitz::Array<int, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_up_ = 0;
    sunpathsfine_up_.reference( blitz::Array<double, 4>(maxlayers_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_up_ = 0;
    ntraversefine_up_.reference( blitz::Array<int, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_up_ = 0;
    xfine_dn_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_dn_ = 0;
    wfine_dn_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_dn_ = 0;
    sunpaths_dn_.reference( blitz::Array<double, 3>(maxlayers_-0+1, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_dn_ = 0;
    ntraverse_dn_.reference( blitz::Array<int, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_dn_ = 0;
    sunpathsfine_dn_.reference( blitz::Array<double, 4>(maxlayers_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_dn_ = 0;
    ntraversefine_dn_.reference( blitz::Array<int, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_dn_ = 0;
    xfine_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_up_p_ = 0;
    wfine_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_up_p_ = 0;
    sunpaths_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_up_p_ = 0;
    ntraverse_up_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_up_p_ = 0;
    sunpathsfine_up_p_.reference( blitz::Array<double, 4>(maxpartials_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_up_p_ = 0;
    ntraversefine_up_p_.reference( blitz::Array<int, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_up_p_ = 0;
    xfine_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_dn_p_ = 0;
    wfine_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_dn_p_ = 0;
    sunpaths_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_dn_p_ = 0;
    ntraverse_dn_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_dn_p_ = 0;
    sunpathsfine_dn_p_.reference( blitz::Array<double, 4>(maxpartials_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_dn_p_ = 0;
    ntraversefine_dn_p_.reference( blitz::Array<int, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_dn_p_ = 0;
    // Initialize type pointers
    
  }

  virtual ~Fo_Scalarss_Rtcalcs_I() = default;

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxpartials() const {
    return maxpartials_;
  }

  

  const int& maxfine() const {
    return maxfine_;
  }

  

  const int& maxmoments_input() const {
    return maxmoments_input_;
  }

  

  const int& max_user_levels() const {
    return max_user_levels_;
  }

  

  const bool& do_deltam_scaling() const {
    return do_deltam_scaling_;
  }

  void do_deltam_scaling(const bool& do_deltam_scaling_in) {
    do_deltam_scaling_ = do_deltam_scaling_in;
  }

  

  const bool& do_phasfunc() const {
    return do_phasfunc_;
  }

  void do_phasfunc(const bool& do_phasfunc_in) {
    do_phasfunc_ = do_phasfunc_in;
  }

  

  const bool& do_partials() const {
    return do_partials_;
  }

  void do_partials(const bool& do_partials_in) {
    do_partials_ = do_partials_in;
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

  

  const double& flux() const {
    return flux_;
  }

  void flux(const double& flux_in) {
    flux_ = flux_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_dn() const {
    return do_sources_dn_;
  }

  void do_sources_dn(const blitz::Array<bool, 2>& do_sources_dn_in) {
    do_sources_dn_ = do_sources_dn_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_dn_p() const {
    return do_sources_dn_p_;
  }

  void do_sources_dn_p(const blitz::Array<bool, 2>& do_sources_dn_p_in) {
    do_sources_dn_p_ = do_sources_dn_p_in;
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

  

  const int& nmoments_input() const {
    return nmoments_input_;
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

  

  const int& npartials() const {
    return npartials_;
  }

  

  const blitz::Array<int, 2>& nfinedivs_p() const {
    return nfinedivs_p_;
  }

  void nfinedivs_p(const blitz::Array<int, 2>& nfinedivs_p_in) {
    nfinedivs_p_ = nfinedivs_p_in;
  }

  

  const blitz::Array<int, 1>& partial_outindex() const {
    return partial_outindex_;
  }

  void partial_outindex(const blitz::Array<int, 1>& partial_outindex_in) {
    partial_outindex_ = partial_outindex_in;
  }

  

  const blitz::Array<bool, 1>& partial_outflag() const {
    return partial_outflag_;
  }

  void partial_outflag(const blitz::Array<bool, 1>& partial_outflag_in) {
    partial_outflag_ = partial_outflag_in;
  }

  

  const blitz::Array<int, 1>& partial_layeridx() const {
    return partial_layeridx_;
  }

  void partial_layeridx(const blitz::Array<int, 1>& partial_layeridx_in) {
    partial_layeridx_ = partial_layeridx_in;
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

  

  const blitz::Array<double, 2>& phasmoms() const {
    return phasmoms_;
  }

  void phasmoms(const blitz::Array<double, 2>& phasmoms_in) {
    phasmoms_ = phasmoms_in;
  }

  

  const blitz::Array<double, 2>& phasfunc_dn() const {
    return phasfunc_dn_;
  }

  void phasfunc_dn(const blitz::Array<double, 2>& phasfunc_dn_in) {
    phasfunc_dn_ = phasfunc_dn_in;
  }

  

  const blitz::Array<double, 1>& mu1() const {
    return mu1_;
  }

  void mu1(const blitz::Array<double, 1>& mu1_in) {
    mu1_ = mu1_in;
  }

  

  const blitz::Array<double, 2>& legpoly_dn() const {
    return legpoly_dn_;
  }

  void legpoly_dn(const blitz::Array<double, 2>& legpoly_dn_in) {
    legpoly_dn_ = legpoly_dn_in;
  }

  

  const blitz::Array<double, 2>& losw_paths() const {
    return losw_paths_;
  }

  void losw_paths(const blitz::Array<double, 2>& losw_paths_in) {
    losw_paths_ = losw_paths_in;
  }

  

  const blitz::Array<double, 2>& losp_paths() const {
    return losp_paths_;
  }

  void losp_paths(const blitz::Array<double, 2>& losp_paths_in) {
    losp_paths_ = losp_paths_in;
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

  

  const blitz::Array<double, 3>& xfine_p() const {
    return xfine_p_;
  }

  void xfine_p(const blitz::Array<double, 3>& xfine_p_in) {
    xfine_p_ = xfine_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_p() const {
    return wfine_p_;
  }

  void wfine_p(const blitz::Array<double, 3>& wfine_p_in) {
    wfine_p_ = wfine_p_in;
  }

  

  const blitz::Array<double, 3>& sunpaths_p() const {
    return sunpaths_p_;
  }

  void sunpaths_p(const blitz::Array<double, 3>& sunpaths_p_in) {
    sunpaths_p_ = sunpaths_p_in;
  }

  

  const blitz::Array<int, 2>& ntraverse_p() const {
    return ntraverse_p_;
  }

  void ntraverse_p(const blitz::Array<int, 2>& ntraverse_p_in) {
    ntraverse_p_ = ntraverse_p_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine_p() const {
    return sunpathsfine_p_;
  }

  void sunpathsfine_p(const blitz::Array<double, 4>& sunpathsfine_p_in) {
    sunpathsfine_p_ = sunpathsfine_p_in;
  }

  

  const blitz::Array<int, 3>& ntraversefine_p() const {
    return ntraversefine_p_;
  }

  void ntraversefine_p(const blitz::Array<int, 3>& ntraversefine_p_in) {
    ntraversefine_p_ = ntraversefine_p_in;
  }

  

  const blitz::Array<double, 2>& intensity_dn() const {
    return intensity_dn_;
  }

  

  const blitz::Array<double, 2>& cumsource_dn() const {
    return cumsource_dn_;
  }

  

  const blitz::Array<double, 2>& lostrans_dn() const {
    return lostrans_dn_;
  }

  

  const bool& do_surface_leaving() const {
    return do_surface_leaving_;
  }

  void do_surface_leaving(const bool& do_surface_leaving_in) {
    do_surface_leaving_ = do_surface_leaving_in;
  }

  

  const bool& do_water_leaving() const {
    return do_water_leaving_;
  }

  void do_water_leaving(const bool& do_water_leaving_in) {
    do_water_leaving_ = do_water_leaving_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_up() const {
    return do_sources_up_;
  }

  void do_sources_up(const blitz::Array<bool, 2>& do_sources_up_in) {
    do_sources_up_ = do_sources_up_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_up_p() const {
    return do_sources_up_p_;
  }

  void do_sources_up_p(const blitz::Array<bool, 2>& do_sources_up_p_in) {
    do_sources_up_p_ = do_sources_up_p_in;
  }

  

  const blitz::Array<double, 2>& phasfunc_up() const {
    return phasfunc_up_;
  }

  void phasfunc_up(const blitz::Array<double, 2>& phasfunc_up_in) {
    phasfunc_up_ = phasfunc_up_in;
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

  

  const blitz::Array<double, 1>& mu0() const {
    return mu0_;
  }

  void mu0(const blitz::Array<double, 1>& mu0_in) {
    mu0_ = mu0_in;
  }

  

  const blitz::Array<double, 2>& legpoly_up() const {
    return legpoly_up_;
  }

  void legpoly_up(const blitz::Array<double, 2>& legpoly_up_in) {
    legpoly_up_ = legpoly_up_in;
  }

  

  const blitz::Array<double, 2>& intensity_up() const {
    return intensity_up_;
  }

  

  const blitz::Array<double, 2>& intensity_db() const {
    return intensity_db_;
  }

  

  const blitz::Array<double, 2>& cumsource_up() const {
    return cumsource_up_;
  }

  

  const blitz::Array<double, 2>& cumtrans() const {
    return cumtrans_;
  }

  

  const blitz::Array<double, 2>& lostrans_up() const {
    return lostrans_up_;
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

  

  const blitz::Array<double, 3>& xfine_up() const {
    return xfine_up_;
  }

  void xfine_up(const blitz::Array<double, 3>& xfine_up_in) {
    xfine_up_ = xfine_up_in;
  }

  

  const blitz::Array<double, 3>& wfine_up() const {
    return wfine_up_;
  }

  void wfine_up(const blitz::Array<double, 3>& wfine_up_in) {
    wfine_up_ = wfine_up_in;
  }

  

  const blitz::Array<double, 3>& sunpaths_up() const {
    return sunpaths_up_;
  }

  void sunpaths_up(const blitz::Array<double, 3>& sunpaths_up_in) {
    sunpaths_up_ = sunpaths_up_in;
  }

  

  const blitz::Array<int, 2>& ntraverse_up() const {
    return ntraverse_up_;
  }

  void ntraverse_up(const blitz::Array<int, 2>& ntraverse_up_in) {
    ntraverse_up_ = ntraverse_up_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine_up() const {
    return sunpathsfine_up_;
  }

  void sunpathsfine_up(const blitz::Array<double, 4>& sunpathsfine_up_in) {
    sunpathsfine_up_ = sunpathsfine_up_in;
  }

  

  const blitz::Array<int, 3>& ntraversefine_up() const {
    return ntraversefine_up_;
  }

  void ntraversefine_up(const blitz::Array<int, 3>& ntraversefine_up_in) {
    ntraversefine_up_ = ntraversefine_up_in;
  }

  

  const blitz::Array<double, 3>& xfine_dn() const {
    return xfine_dn_;
  }

  void xfine_dn(const blitz::Array<double, 3>& xfine_dn_in) {
    xfine_dn_ = xfine_dn_in;
  }

  

  const blitz::Array<double, 3>& wfine_dn() const {
    return wfine_dn_;
  }

  void wfine_dn(const blitz::Array<double, 3>& wfine_dn_in) {
    wfine_dn_ = wfine_dn_in;
  }

  

  const blitz::Array<double, 3>& sunpaths_dn() const {
    return sunpaths_dn_;
  }

  void sunpaths_dn(const blitz::Array<double, 3>& sunpaths_dn_in) {
    sunpaths_dn_ = sunpaths_dn_in;
  }

  

  const blitz::Array<int, 2>& ntraverse_dn() const {
    return ntraverse_dn_;
  }

  void ntraverse_dn(const blitz::Array<int, 2>& ntraverse_dn_in) {
    ntraverse_dn_ = ntraverse_dn_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine_dn() const {
    return sunpathsfine_dn_;
  }

  void sunpathsfine_dn(const blitz::Array<double, 4>& sunpathsfine_dn_in) {
    sunpathsfine_dn_ = sunpathsfine_dn_in;
  }

  

  const blitz::Array<int, 3>& ntraversefine_dn() const {
    return ntraversefine_dn_;
  }

  void ntraversefine_dn(const blitz::Array<int, 3>& ntraversefine_dn_in) {
    ntraversefine_dn_ = ntraversefine_dn_in;
  }

  

  const blitz::Array<double, 3>& xfine_up_p() const {
    return xfine_up_p_;
  }

  void xfine_up_p(const blitz::Array<double, 3>& xfine_up_p_in) {
    xfine_up_p_ = xfine_up_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_up_p() const {
    return wfine_up_p_;
  }

  void wfine_up_p(const blitz::Array<double, 3>& wfine_up_p_in) {
    wfine_up_p_ = wfine_up_p_in;
  }

  

  const blitz::Array<double, 3>& sunpaths_up_p() const {
    return sunpaths_up_p_;
  }

  void sunpaths_up_p(const blitz::Array<double, 3>& sunpaths_up_p_in) {
    sunpaths_up_p_ = sunpaths_up_p_in;
  }

  

  const blitz::Array<int, 2>& ntraverse_up_p() const {
    return ntraverse_up_p_;
  }

  void ntraverse_up_p(const blitz::Array<int, 2>& ntraverse_up_p_in) {
    ntraverse_up_p_ = ntraverse_up_p_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine_up_p() const {
    return sunpathsfine_up_p_;
  }

  void sunpathsfine_up_p(const blitz::Array<double, 4>& sunpathsfine_up_p_in) {
    sunpathsfine_up_p_ = sunpathsfine_up_p_in;
  }

  

  const blitz::Array<int, 3>& ntraversefine_up_p() const {
    return ntraversefine_up_p_;
  }

  void ntraversefine_up_p(const blitz::Array<int, 3>& ntraversefine_up_p_in) {
    ntraversefine_up_p_ = ntraversefine_up_p_in;
  }

  

  const blitz::Array<double, 3>& xfine_dn_p() const {
    return xfine_dn_p_;
  }

  void xfine_dn_p(const blitz::Array<double, 3>& xfine_dn_p_in) {
    xfine_dn_p_ = xfine_dn_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_dn_p() const {
    return wfine_dn_p_;
  }

  void wfine_dn_p(const blitz::Array<double, 3>& wfine_dn_p_in) {
    wfine_dn_p_ = wfine_dn_p_in;
  }

  

  const blitz::Array<double, 3>& sunpaths_dn_p() const {
    return sunpaths_dn_p_;
  }

  void sunpaths_dn_p(const blitz::Array<double, 3>& sunpaths_dn_p_in) {
    sunpaths_dn_p_ = sunpaths_dn_p_in;
  }

  

  const blitz::Array<int, 2>& ntraverse_dn_p() const {
    return ntraverse_dn_p_;
  }

  void ntraverse_dn_p(const blitz::Array<int, 2>& ntraverse_dn_p_in) {
    ntraverse_dn_p_ = ntraverse_dn_p_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine_dn_p() const {
    return sunpathsfine_dn_p_;
  }

  void sunpathsfine_dn_p(const blitz::Array<double, 4>& sunpathsfine_dn_p_in) {
    sunpathsfine_dn_p_ = sunpathsfine_dn_p_in;
  }

  

  const blitz::Array<int, 3>& ntraversefine_dn_p() const {
    return ntraversefine_dn_p_;
  }

  void ntraversefine_dn_p(const blitz::Array<int, 3>& ntraversefine_dn_p_in) {
    ntraversefine_dn_p_ = ntraversefine_dn_p_in;
  }

  

  
  void ss_integral_i_dn() {
    
    
    fo_scalarss_rtcalcs_i_m_ss_integral_i_dn_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &maxmoments_input_, &max_user_levels_, &do_deltam_scaling_, &do_phasfunc_, &do_partials_, &do_planpar_, &do_enhanced_ps_, &flux_, do_sources_dn_.dataFirst(), do_sources_dn_p_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &nmoments_input_, &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), phasmoms_.dataFirst(), phasfunc_dn_.dataFirst(), mu1_.dataFirst(), legpoly_dn_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), sunpaths_.dataFirst(), ntraverse_.dataFirst(), sunpathsfine_.dataFirst(), ntraversefine_.dataFirst(), xfine_p_.dataFirst(), wfine_p_.dataFirst(), sunpaths_p_.dataFirst(), ntraverse_p_.dataFirst(), sunpathsfine_p_.dataFirst(), ntraversefine_p_.dataFirst(), intensity_dn_.dataFirst(), cumsource_dn_.dataFirst(), lostrans_dn_.dataFirst());
    
  }
void ss_integral_i_up() {
    
    
    fo_scalarss_rtcalcs_i_m_ss_integral_i_up_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &maxmoments_input_, &max_user_levels_, &do_deltam_scaling_, &do_phasfunc_, &do_surface_leaving_, &do_water_leaving_, &do_partials_, &do_planpar_, &do_enhanced_ps_, &flux_, do_sources_up_.dataFirst(), do_sources_up_p_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &nmoments_input_, &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), phasmoms_.dataFirst(), phasfunc_up_.dataFirst(), reflec_.dataFirst(), slterm_.dataFirst(), mu0_.dataFirst(), mu1_.dataFirst(), legpoly_up_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), sunpaths_.dataFirst(), ntraverse_.dataFirst(), sunpathsfine_.dataFirst(), ntraversefine_.dataFirst(), xfine_p_.dataFirst(), wfine_p_.dataFirst(), sunpaths_p_.dataFirst(), ntraverse_p_.dataFirst(), sunpathsfine_p_.dataFirst(), ntraversefine_p_.dataFirst(), intensity_up_.dataFirst(), intensity_db_.dataFirst(), cumsource_up_.dataFirst(), cumtrans_.dataFirst(), lostrans_up_.dataFirst());
    
  }
void ss_integral_i_updn() {
    
    
    fo_scalarss_rtcalcs_i_m_ss_integral_i_updn_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &maxmoments_input_, &max_user_levels_, &do_upwelling_, &do_dnwelling_, &do_deltam_scaling_, &do_phasfunc_, &do_surface_leaving_, &do_water_leaving_, &do_partials_, &do_planpar_, &do_enhanced_ps_, do_sources_up_.dataFirst(), do_sources_up_p_.dataFirst(), do_sources_dn_.dataFirst(), do_sources_dn_p_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &nmoments_input_, &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), &flux_, extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), phasmoms_.dataFirst(), phasfunc_up_.dataFirst(), phasfunc_dn_.dataFirst(), reflec_.dataFirst(), slterm_.dataFirst(), mu0_.dataFirst(), mu1_.dataFirst(), legpoly_up_.dataFirst(), legpoly_dn_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_up_.dataFirst(), wfine_up_.dataFirst(), sunpaths_up_.dataFirst(), ntraverse_up_.dataFirst(), sunpathsfine_up_.dataFirst(), ntraversefine_up_.dataFirst(), xfine_dn_.dataFirst(), wfine_dn_.dataFirst(), sunpaths_dn_.dataFirst(), ntraverse_dn_.dataFirst(), sunpathsfine_dn_.dataFirst(), ntraversefine_dn_.dataFirst(), xfine_up_p_.dataFirst(), wfine_up_p_.dataFirst(), sunpaths_up_p_.dataFirst(), ntraverse_up_p_.dataFirst(), sunpathsfine_up_p_.dataFirst(), ntraversefine_up_p_.dataFirst(), xfine_dn_p_.dataFirst(), wfine_dn_p_.dataFirst(), sunpaths_dn_p_.dataFirst(), ntraverse_dn_p_.dataFirst(), sunpathsfine_dn_p_.dataFirst(), ntraversefine_dn_p_.dataFirst(), intensity_up_.dataFirst(), intensity_db_.dataFirst(), cumsource_up_.dataFirst(), cumtrans_.dataFirst(), lostrans_up_.dataFirst(), intensity_dn_.dataFirst(), cumsource_dn_.dataFirst(), lostrans_dn_.dataFirst());
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "Fo_Scalarss_Rtcalcs_I:" << std::endl
      << "          maxgeoms: " << maxgeoms()  << std::endl
      << "         maxlayers: " << maxlayers()  << std::endl
      << "       maxpartials: " << maxpartials()  << std::endl
      << "           maxfine: " << maxfine()  << std::endl
      << "  maxmoments_input: " << maxmoments_input()  << std::endl
      << "   max_user_levels: " << max_user_levels()  << std::endl
      << " do_deltam_scaling: " << do_deltam_scaling()  << std::endl
      << "       do_phasfunc: " << do_phasfunc()  << std::endl
      << "       do_partials: " << do_partials()  << std::endl
      << "        do_planpar: " << do_planpar()  << std::endl
      << "    do_enhanced_ps: " << do_enhanced_ps()  << std::endl
      << "              flux: " << flux()  << std::endl
      << "     do_sources_dn: " << std::endl << do_sources_dn()  << std::endl
      << "   do_sources_dn_p: " << std::endl << do_sources_dn_p()  << std::endl
      << "            ngeoms: " << ngeoms()  << std::endl
      << "           nlayers: " << nlayers()  << std::endl
      << "         nfinedivs: " << std::endl << nfinedivs()  << std::endl
      << "    nmoments_input: " << nmoments_input()  << std::endl
      << "     n_user_levels: " << n_user_levels()  << std::endl
      << "       user_levels: " << std::endl << user_levels()  << std::endl
      << "         npartials: " << npartials()  << std::endl
      << "       nfinedivs_p: " << std::endl << nfinedivs_p()  << std::endl
      << "  partial_outindex: " << std::endl << partial_outindex()  << std::endl
      << "   partial_outflag: " << std::endl << partial_outflag()  << std::endl
      << "  partial_layeridx: " << std::endl << partial_layeridx()  << std::endl
      << "        extinction: " << std::endl << extinction()  << std::endl
      << "           deltaus: " << std::endl << deltaus()  << std::endl
      << "             omega: " << std::endl << omega()  << std::endl
      << "          truncfac: " << std::endl << truncfac()  << std::endl
      << "          phasmoms: " << std::endl << phasmoms()  << std::endl
      << "       phasfunc_dn: " << std::endl << phasfunc_dn()  << std::endl
      << "               mu1: " << std::endl << mu1()  << std::endl
      << "        legpoly_dn: " << std::endl << legpoly_dn()  << std::endl
      << "        losw_paths: " << std::endl << losw_paths()  << std::endl
      << "        losp_paths: " << std::endl << losp_paths()  << std::endl
      << "             xfine: " << std::endl << xfine()  << std::endl
      << "             wfine: " << std::endl << wfine()  << std::endl
      << "          sunpaths: " << std::endl << sunpaths()  << std::endl
      << "         ntraverse: " << std::endl << ntraverse()  << std::endl
      << "      sunpathsfine: " << std::endl << sunpathsfine()  << std::endl
      << "     ntraversefine: " << std::endl << ntraversefine()  << std::endl
      << "           xfine_p: " << std::endl << xfine_p()  << std::endl
      << "           wfine_p: " << std::endl << wfine_p()  << std::endl
      << "        sunpaths_p: " << std::endl << sunpaths_p()  << std::endl
      << "       ntraverse_p: " << std::endl << ntraverse_p()  << std::endl
      << "    sunpathsfine_p: " << std::endl << sunpathsfine_p()  << std::endl
      << "   ntraversefine_p: " << std::endl << ntraversefine_p()  << std::endl
      << "      intensity_dn: " << std::endl << intensity_dn()  << std::endl
      << "      cumsource_dn: " << std::endl << cumsource_dn()  << std::endl
      << "       lostrans_dn: " << std::endl << lostrans_dn()  << std::endl
      << "do_surface_leaving: " << do_surface_leaving()  << std::endl
      << "  do_water_leaving: " << do_water_leaving()  << std::endl
      << "     do_sources_up: " << std::endl << do_sources_up()  << std::endl
      << "   do_sources_up_p: " << std::endl << do_sources_up_p()  << std::endl
      << "       phasfunc_up: " << std::endl << phasfunc_up()  << std::endl
      << "            reflec: " << std::endl << reflec()  << std::endl
      << "            slterm: " << std::endl << slterm()  << std::endl
      << "               mu0: " << std::endl << mu0()  << std::endl
      << "        legpoly_up: " << std::endl << legpoly_up()  << std::endl
      << "      intensity_up: " << std::endl << intensity_up()  << std::endl
      << "      intensity_db: " << std::endl << intensity_db()  << std::endl
      << "      cumsource_up: " << std::endl << cumsource_up()  << std::endl
      << "          cumtrans: " << std::endl << cumtrans()  << std::endl
      << "       lostrans_up: " << std::endl << lostrans_up()  << std::endl
      << "      do_upwelling: " << do_upwelling()  << std::endl
      << "      do_dnwelling: " << do_dnwelling()  << std::endl
      << "          xfine_up: " << std::endl << xfine_up()  << std::endl
      << "          wfine_up: " << std::endl << wfine_up()  << std::endl
      << "       sunpaths_up: " << std::endl << sunpaths_up()  << std::endl
      << "      ntraverse_up: " << std::endl << ntraverse_up()  << std::endl
      << "   sunpathsfine_up: " << std::endl << sunpathsfine_up()  << std::endl
      << "  ntraversefine_up: " << std::endl << ntraversefine_up()  << std::endl
      << "          xfine_dn: " << std::endl << xfine_dn()  << std::endl
      << "          wfine_dn: " << std::endl << wfine_dn()  << std::endl
      << "       sunpaths_dn: " << std::endl << sunpaths_dn()  << std::endl
      << "      ntraverse_dn: " << std::endl << ntraverse_dn()  << std::endl
      << "   sunpathsfine_dn: " << std::endl << sunpathsfine_dn()  << std::endl
      << "  ntraversefine_dn: " << std::endl << ntraversefine_dn()  << std::endl
      << "        xfine_up_p: " << std::endl << xfine_up_p()  << std::endl
      << "        wfine_up_p: " << std::endl << wfine_up_p()  << std::endl
      << "     sunpaths_up_p: " << std::endl << sunpaths_up_p()  << std::endl
      << "    ntraverse_up_p: " << std::endl << ntraverse_up_p()  << std::endl
      << " sunpathsfine_up_p: " << std::endl << sunpathsfine_up_p()  << std::endl
      << "ntraversefine_up_p: " << std::endl << ntraversefine_up_p()  << std::endl
      << "        xfine_dn_p: " << std::endl << xfine_dn_p()  << std::endl
      << "        wfine_dn_p: " << std::endl << wfine_dn_p()  << std::endl
      << "     sunpaths_dn_p: " << std::endl << sunpaths_dn_p()  << std::endl
      << "    ntraverse_dn_p: " << std::endl << ntraverse_dn_p()  << std::endl
      << " sunpathsfine_dn_p: " << std::endl << sunpathsfine_dn_p()  << std::endl
      << "ntraversefine_dn_p: " << std::endl << ntraversefine_dn_p()  << std::endl;

  }

private:
  int maxgeoms_;
  int maxlayers_;
  int maxpartials_;
  int maxfine_;
  int maxmoments_input_;
  int max_user_levels_;
  bool do_deltam_scaling_;
  bool do_phasfunc_;
  bool do_partials_;
  bool do_planpar_;
  bool do_enhanced_ps_;
  double flux_;
  blitz::Array<bool, 2> do_sources_dn_;
  blitz::Array<bool, 2> do_sources_dn_p_;
  int ngeoms_;
  int nlayers_;
  blitz::Array<int, 2> nfinedivs_;
  int nmoments_input_;
  int n_user_levels_;
  blitz::Array<int, 1> user_levels_;
  int npartials_;
  blitz::Array<int, 2> nfinedivs_p_;
  blitz::Array<int, 1> partial_outindex_;
  blitz::Array<bool, 1> partial_outflag_;
  blitz::Array<int, 1> partial_layeridx_;
  blitz::Array<double, 1> extinction_;
  blitz::Array<double, 1> deltaus_;
  blitz::Array<double, 1> omega_;
  blitz::Array<double, 1> truncfac_;
  blitz::Array<double, 2> phasmoms_;
  blitz::Array<double, 2> phasfunc_dn_;
  blitz::Array<double, 1> mu1_;
  blitz::Array<double, 2> legpoly_dn_;
  blitz::Array<double, 2> losw_paths_;
  blitz::Array<double, 2> losp_paths_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> sunpaths_;
  blitz::Array<int, 2> ntraverse_;
  blitz::Array<double, 4> sunpathsfine_;
  blitz::Array<int, 3> ntraversefine_;
  blitz::Array<double, 3> xfine_p_;
  blitz::Array<double, 3> wfine_p_;
  blitz::Array<double, 3> sunpaths_p_;
  blitz::Array<int, 2> ntraverse_p_;
  blitz::Array<double, 4> sunpathsfine_p_;
  blitz::Array<int, 3> ntraversefine_p_;
  blitz::Array<double, 2> intensity_dn_;
  blitz::Array<double, 2> cumsource_dn_;
  blitz::Array<double, 2> lostrans_dn_;
  bool do_surface_leaving_;
  bool do_water_leaving_;
  blitz::Array<bool, 2> do_sources_up_;
  blitz::Array<bool, 2> do_sources_up_p_;
  blitz::Array<double, 2> phasfunc_up_;
  blitz::Array<double, 1> reflec_;
  blitz::Array<double, 1> slterm_;
  blitz::Array<double, 1> mu0_;
  blitz::Array<double, 2> legpoly_up_;
  blitz::Array<double, 2> intensity_up_;
  blitz::Array<double, 2> intensity_db_;
  blitz::Array<double, 2> cumsource_up_;
  blitz::Array<double, 2> cumtrans_;
  blitz::Array<double, 2> lostrans_up_;
  bool do_upwelling_;
  bool do_dnwelling_;
  blitz::Array<double, 3> xfine_up_;
  blitz::Array<double, 3> wfine_up_;
  blitz::Array<double, 3> sunpaths_up_;
  blitz::Array<int, 2> ntraverse_up_;
  blitz::Array<double, 4> sunpathsfine_up_;
  blitz::Array<int, 3> ntraversefine_up_;
  blitz::Array<double, 3> xfine_dn_;
  blitz::Array<double, 3> wfine_dn_;
  blitz::Array<double, 3> sunpaths_dn_;
  blitz::Array<int, 2> ntraverse_dn_;
  blitz::Array<double, 4> sunpathsfine_dn_;
  blitz::Array<int, 3> ntraversefine_dn_;
  blitz::Array<double, 3> xfine_up_p_;
  blitz::Array<double, 3> wfine_up_p_;
  blitz::Array<double, 3> sunpaths_up_p_;
  blitz::Array<int, 2> ntraverse_up_p_;
  blitz::Array<double, 4> sunpathsfine_up_p_;
  blitz::Array<int, 3> ntraversefine_up_p_;
  blitz::Array<double, 3> xfine_dn_p_;
  blitz::Array<double, 3> wfine_dn_p_;
  blitz::Array<double, 3> sunpaths_dn_p_;
  blitz::Array<int, 2> ntraverse_dn_p_;
  blitz::Array<double, 4> sunpathsfine_dn_p_;
  blitz::Array<int, 3> ntraversefine_dn_p_;

  // Serialization support
  Fo_Scalarss_Rtcalcs_I() {}
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "fo_scalarss_rtcalcs_ilps_m" in file: "FO_ScalarSS_RTCalcs_ILPS.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_dn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* maxmoments_input_in, const int* max_user_levels_in, const int* max_atmoswfs_in, const bool* do_deltam_scaling_in, const bool* do_phasfunc_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const double* flux_in, const bool* do_sources_dn_in, const bool* do_sources_dn_p_in, const bool* do_profilewfs_in, const bool* lvaryflags_in, const int* lvarynums_in, const bool* lvarymoms_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* nmoments_input_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* phasmoms_in, const double* phasfunc_dn_in, const double* l_extinction_in, const double* l_deltaus_in, const double* l_omega_in, const double* l_truncfac_in, const double* l_phasmoms_in, const double* l_phasfunc_dn_in, const double* mu1_in, const double* legpoly_dn_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_in, const double* wfine_in, const double* sunpaths_in, const int* ntraverse_in, const double* sunpathsfine_in, const int* ntraversefine_in, const double* xfine_p_in, const double* wfine_p_in, const double* sunpaths_p_in, const int* ntraverse_p_in, const double* sunpathsfine_p_in, const int* ntraversefine_p_in, const double* intensity_dn_in, const double* lp_jacobians_dn_in, const double* lostrans_dn_in, const double* lp_lostrans_dn_in);
  void fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_up_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* maxmoments_input_in, const int* max_user_levels_in, const int* max_atmoswfs_in, const int* max_surfacewfs_in, const int* max_sleavewfs_in, const bool* do_deltam_scaling_in, const bool* do_phasfunc_in, const bool* do_surface_leaving_in, const bool* do_water_leaving_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const double* flux_in, const bool* do_sources_up_in, const bool* do_sources_up_p_in, const bool* do_profilewfs_in, const bool* do_surfacewfs_in, const bool* do_sleavewfs_in, const int* n_reflecwfs_in, const int* n_sleavewfs_in, const int* n_surfacewfs_in, const bool* lvaryflags_in, const int* lvarynums_in, const bool* lvarymoms_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* nmoments_input_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* phasmoms_in, const double* phasfunc_up_in, const double* reflec_in, const double* slterm_in, const double* l_extinction_in, const double* l_deltaus_in, const double* l_omega_in, const double* l_truncfac_in, const double* l_phasmoms_in, const double* l_phasfunc_up_in, const double* ls_reflec_in, const double* lssl_slterm_in, const double* mu0_in, const double* mu1_in, const double* legpoly_up_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_in, const double* wfine_in, const double* sunpaths_in, const int* ntraverse_in, const double* sunpathsfine_in, const int* ntraversefine_in, const double* xfine_p_in, const double* wfine_p_in, const double* sunpaths_p_in, const int* ntraverse_p_in, const double* sunpathsfine_p_in, const int* ntraversefine_p_in, const double* intensity_up_in, const double* intensity_db_in, const double* lp_jacobians_up_in, const double* lp_jacobians_db_in, const double* ls_jacobians_db_in, const double* cumtrans_in, const double* lostrans_up_in, const double* lp_cumtrans_in, const double* lp_lostrans_up_in);
  void fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_updn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* maxmoments_input_in, const int* max_user_levels_in, const int* max_atmoswfs_in, const int* max_surfacewfs_in, const int* max_sleavewfs_in, const bool* do_upwelling_in, const bool* do_dnwelling_in, const bool* do_deltam_scaling_in, const bool* do_phasfunc_in, const bool* do_surface_leaving_in, const bool* do_water_leaving_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const bool* do_sources_up_in, const bool* do_sources_up_p_in, const bool* do_sources_dn_in, const bool* do_sources_dn_p_in, const bool* do_profilewfs_in, const bool* do_surfacewfs_in, const bool* do_sleavewfs_in, const int* n_reflecwfs_in, const int* n_sleavewfs_in, const int* n_surfacewfs_in, const bool* lvaryflags_in, const int* lvarynums_in, const bool* lvarymoms_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* nmoments_input_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* flux_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* phasmoms_in, const double* phasfunc_up_in, const double* phasfunc_dn_in, const double* reflec_in, const double* slterm_in, const double* ls_reflec_in, const double* lssl_slterm_in, const double* l_extinction_in, const double* l_deltaus_in, const double* l_omega_in, const double* l_truncfac_in, const double* l_phasmoms_in, const double* l_phasfunc_up_in, const double* l_phasfunc_dn_in, const double* mu0_in, const double* mu1_in, const double* legpoly_up_in, const double* legpoly_dn_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_up_in, const double* wfine_up_in, const double* sunpaths_up_in, const int* ntraverse_up_in, const double* sunpathsfine_up_in, const int* ntraversefine_up_in, const double* xfine_dn_in, const double* wfine_dn_in, const double* sunpaths_dn_in, const int* ntraverse_dn_in, const double* sunpathsfine_dn_in, const int* ntraversefine_dn_in, const double* xfine_up_p_in, const double* wfine_up_p_in, const double* sunpaths_up_p_in, const int* ntraverse_up_p_in, const double* sunpathsfine_up_p_in, const int* ntraversefine_up_p_in, const double* xfine_dn_p_in, const double* wfine_dn_p_in, const double* sunpaths_dn_p_in, const int* ntraverse_dn_p_in, const double* sunpathsfine_dn_p_in, const int* ntraversefine_dn_p_in, const double* intensity_up_in, const double* intensity_db_in, const double* lp_jacobians_up_in, const double* lp_jacobians_db_in, const double* ls_jacobians_db_in, const double* lostrans_up_in, const double* lp_lostrans_up_in, const double* intensity_dn_in, const double* lp_jacobians_dn_in, const double* cumtrans_in, const double* lp_cumtrans_in, const double* lostrans_dn_in, const double* lp_lostrans_dn_in);
}

class Fo_Scalarss_Rtcalcs_Ilps : public Printable<Fo_Scalarss_Rtcalcs_Ilps> {

public:
  Fo_Scalarss_Rtcalcs_Ilps(const int& maxgeoms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& maxmoments_input_in, const int& max_user_levels_in, const int& max_atmoswfs_in, const int& ngeoms_in, const int& nlayers_in, const int& nmoments_input_in, const int& n_user_levels_in, const int& npartials_in, const int& max_surfacewfs_in, const int& max_sleavewfs_in, const int& n_sleavewfs_in, const int& n_surfacewfs_in) : maxgeoms_(maxgeoms_in), maxlayers_(maxlayers_in), maxpartials_(maxpartials_in), maxfine_(maxfine_in), maxmoments_input_(maxmoments_input_in), max_user_levels_(max_user_levels_in), max_atmoswfs_(max_atmoswfs_in), ngeoms_(ngeoms_in), nlayers_(nlayers_in), nmoments_input_(nmoments_input_in), n_user_levels_(n_user_levels_in), npartials_(npartials_in), max_surfacewfs_(max_surfacewfs_in), max_sleavewfs_(max_sleavewfs_in), n_sleavewfs_(n_sleavewfs_in), n_surfacewfs_(n_surfacewfs_in) 
  { 
    do_deltam_scaling_ = false;
    do_phasfunc_ = false;
    do_partials_ = false;
    do_planpar_ = false;
    do_enhanced_ps_ = false;
    flux_ = 0;
    do_sources_dn_.reference( blitz::Array<bool, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_dn_ = false;
    do_sources_dn_p_.reference( blitz::Array<bool, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_dn_p_ = false;
    do_profilewfs_ = false;
    lvaryflags_.reference( blitz::Array<bool, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    lvaryflags_ = false;
    lvarynums_.reference( blitz::Array<int, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    lvarynums_ = 0;
    lvarymoms_.reference( blitz::Array<bool, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    lvarymoms_ = false;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    user_levels_.reference( blitz::Array<int, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    user_levels_ = 0;
    nfinedivs_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_p_ = 0;
    partial_outindex_.reference( blitz::Array<int, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    partial_outindex_ = 0;
    partial_outflag_.reference( blitz::Array<bool, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    partial_outflag_ = false;
    partial_layeridx_.reference( blitz::Array<int, 1>(maxpartials_, blitz::ColumnMajorArray<1>()) );
    partial_layeridx_ = 0;
    extinction_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    extinction_ = 0;
    deltaus_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    deltaus_ = 0;
    omega_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    omega_ = 0;
    truncfac_.reference( blitz::Array<double, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    truncfac_ = 0;
    phasmoms_.reference( blitz::Array<double, 2>(maxlayers_, maxmoments_input_-0+1, blitz::ColumnMajorArray<2>()) );
    phasmoms_ = 0;
    phasfunc_dn_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    phasfunc_dn_ = 0;
    l_extinction_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_extinction_ = 0;
    l_deltaus_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_deltaus_ = 0;
    l_omega_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_omega_ = 0;
    l_truncfac_.reference( blitz::Array<double, 2>(maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<2>()) );
    l_truncfac_ = 0;
    l_phasmoms_.reference( blitz::Array<double, 3>(maxlayers_, maxmoments_input_-0+1, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    l_phasmoms_ = 0;
    l_phasfunc_dn_.reference( blitz::Array<double, 3>(maxlayers_, maxgeoms_, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    l_phasfunc_dn_ = 0;
    mu1_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu1_ = 0;
    legpoly_dn_.reference( blitz::Array<double, 2>(maxmoments_input_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    legpoly_dn_ = 0;
    losw_paths_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losw_paths_ = 0;
    losp_paths_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losp_paths_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    sunpaths_.reference( blitz::Array<double, 3>(maxlayers_-0+1, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_ = 0;
    ntraverse_.reference( blitz::Array<int, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_ = 0;
    sunpathsfine_.reference( blitz::Array<double, 4>(maxlayers_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_ = 0;
    ntraversefine_.reference( blitz::Array<int, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_ = 0;
    xfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_p_ = 0;
    wfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_p_ = 0;
    sunpaths_p_.reference( blitz::Array<double, 3>(maxpartials_, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_p_ = 0;
    ntraverse_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_p_ = 0;
    sunpathsfine_p_.reference( blitz::Array<double, 4>(maxpartials_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_p_ = 0;
    ntraversefine_p_.reference( blitz::Array<int, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_p_ = 0;
    intensity_dn_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dn_ = 0;
    lp_jacobians_dn_.reference( blitz::Array<double, 4>(max_user_levels_, maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    lp_jacobians_dn_ = 0;
    lostrans_dn_.reference( blitz::Array<double, 2>(maxgeoms_, maxlayers_, blitz::ColumnMajorArray<2>()) );
    lostrans_dn_ = 0;
    lp_lostrans_dn_.reference( blitz::Array<double, 3>(maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    lp_lostrans_dn_ = 0;
    do_surface_leaving_ = false;
    do_water_leaving_ = false;
    do_sources_up_.reference( blitz::Array<bool, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_up_ = false;
    do_sources_up_p_.reference( blitz::Array<bool, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_up_p_ = false;
    do_surfacewfs_ = false;
    do_sleavewfs_ = false;
    n_reflecwfs_ = 0;
    phasfunc_up_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    phasfunc_up_ = 0;
    reflec_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    reflec_ = 0;
    slterm_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    slterm_ = 0;
    l_phasfunc_up_.reference( blitz::Array<double, 3>(maxlayers_, maxgeoms_, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    l_phasfunc_up_ = 0;
    ls_reflec_.reference( blitz::Array<double, 2>(maxgeoms_, max_surfacewfs_, blitz::ColumnMajorArray<2>()) );
    ls_reflec_ = 0;
    lssl_slterm_.reference( blitz::Array<double, 2>(maxgeoms_, max_sleavewfs_, blitz::ColumnMajorArray<2>()) );
    lssl_slterm_ = 0;
    mu0_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    mu0_ = 0;
    legpoly_up_.reference( blitz::Array<double, 2>(maxmoments_input_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    legpoly_up_ = 0;
    intensity_up_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_up_ = 0;
    intensity_db_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_db_ = 0;
    lp_jacobians_up_.reference( blitz::Array<double, 4>(max_user_levels_, maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    lp_jacobians_up_ = 0;
    lp_jacobians_db_.reference( blitz::Array<double, 4>(max_user_levels_, maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    lp_jacobians_db_ = 0;
    ls_jacobians_db_.reference( blitz::Array<double, 3>(max_user_levels_, maxgeoms_, max_surfacewfs_, blitz::ColumnMajorArray<3>()) );
    ls_jacobians_db_ = 0;
    cumtrans_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cumtrans_ = 0;
    lostrans_up_.reference( blitz::Array<double, 2>(maxgeoms_, maxlayers_, blitz::ColumnMajorArray<2>()) );
    lostrans_up_ = 0;
    lp_cumtrans_.reference( blitz::Array<double, 4>(max_user_levels_, maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    lp_cumtrans_ = 0;
    lp_lostrans_up_.reference( blitz::Array<double, 3>(maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    lp_lostrans_up_ = 0;
    do_upwelling_ = false;
    do_dnwelling_ = false;
    xfine_up_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_up_ = 0;
    wfine_up_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_up_ = 0;
    sunpaths_up_.reference( blitz::Array<double, 3>(maxlayers_-0+1, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_up_ = 0;
    ntraverse_up_.reference( blitz::Array<int, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_up_ = 0;
    sunpathsfine_up_.reference( blitz::Array<double, 4>(maxlayers_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_up_ = 0;
    ntraversefine_up_.reference( blitz::Array<int, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_up_ = 0;
    xfine_dn_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_dn_ = 0;
    wfine_dn_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_dn_ = 0;
    sunpaths_dn_.reference( blitz::Array<double, 3>(maxlayers_-0+1, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_dn_ = 0;
    ntraverse_dn_.reference( blitz::Array<int, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_dn_ = 0;
    sunpathsfine_dn_.reference( blitz::Array<double, 4>(maxlayers_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_dn_ = 0;
    ntraversefine_dn_.reference( blitz::Array<int, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_dn_ = 0;
    xfine_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_up_p_ = 0;
    wfine_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_up_p_ = 0;
    sunpaths_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_up_p_ = 0;
    ntraverse_up_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_up_p_ = 0;
    sunpathsfine_up_p_.reference( blitz::Array<double, 4>(maxpartials_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_up_p_ = 0;
    ntraversefine_up_p_.reference( blitz::Array<int, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_up_p_ = 0;
    xfine_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_dn_p_ = 0;
    wfine_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_dn_p_ = 0;
    sunpaths_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxlayers_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    sunpaths_dn_p_ = 0;
    ntraverse_dn_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ntraverse_dn_p_ = 0;
    sunpathsfine_dn_p_.reference( blitz::Array<double, 4>(maxpartials_, maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<4>()) );
    sunpathsfine_dn_p_ = 0;
    ntraversefine_dn_p_.reference( blitz::Array<int, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    ntraversefine_dn_p_ = 0;
    // Initialize type pointers
    
  }

  virtual ~Fo_Scalarss_Rtcalcs_Ilps() = default;

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxpartials() const {
    return maxpartials_;
  }

  

  const int& maxfine() const {
    return maxfine_;
  }

  

  const int& maxmoments_input() const {
    return maxmoments_input_;
  }

  

  const int& max_user_levels() const {
    return max_user_levels_;
  }

  

  const int& max_atmoswfs() const {
    return max_atmoswfs_;
  }

  

  const bool& do_deltam_scaling() const {
    return do_deltam_scaling_;
  }

  void do_deltam_scaling(const bool& do_deltam_scaling_in) {
    do_deltam_scaling_ = do_deltam_scaling_in;
  }

  

  const bool& do_phasfunc() const {
    return do_phasfunc_;
  }

  void do_phasfunc(const bool& do_phasfunc_in) {
    do_phasfunc_ = do_phasfunc_in;
  }

  

  const bool& do_partials() const {
    return do_partials_;
  }

  void do_partials(const bool& do_partials_in) {
    do_partials_ = do_partials_in;
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

  

  const double& flux() const {
    return flux_;
  }

  void flux(const double& flux_in) {
    flux_ = flux_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_dn() const {
    return do_sources_dn_;
  }

  void do_sources_dn(const blitz::Array<bool, 2>& do_sources_dn_in) {
    do_sources_dn_ = do_sources_dn_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_dn_p() const {
    return do_sources_dn_p_;
  }

  void do_sources_dn_p(const blitz::Array<bool, 2>& do_sources_dn_p_in) {
    do_sources_dn_p_ = do_sources_dn_p_in;
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

  

  const blitz::Array<bool, 2>& lvarymoms() const {
    return lvarymoms_;
  }

  void lvarymoms(const blitz::Array<bool, 2>& lvarymoms_in) {
    lvarymoms_ = lvarymoms_in;
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

  

  const int& nmoments_input() const {
    return nmoments_input_;
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

  

  const int& npartials() const {
    return npartials_;
  }

  

  const blitz::Array<int, 2>& nfinedivs_p() const {
    return nfinedivs_p_;
  }

  void nfinedivs_p(const blitz::Array<int, 2>& nfinedivs_p_in) {
    nfinedivs_p_ = nfinedivs_p_in;
  }

  

  const blitz::Array<int, 1>& partial_outindex() const {
    return partial_outindex_;
  }

  void partial_outindex(const blitz::Array<int, 1>& partial_outindex_in) {
    partial_outindex_ = partial_outindex_in;
  }

  

  const blitz::Array<bool, 1>& partial_outflag() const {
    return partial_outflag_;
  }

  void partial_outflag(const blitz::Array<bool, 1>& partial_outflag_in) {
    partial_outflag_ = partial_outflag_in;
  }

  

  const blitz::Array<int, 1>& partial_layeridx() const {
    return partial_layeridx_;
  }

  void partial_layeridx(const blitz::Array<int, 1>& partial_layeridx_in) {
    partial_layeridx_ = partial_layeridx_in;
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

  

  const blitz::Array<double, 2>& phasmoms() const {
    return phasmoms_;
  }

  void phasmoms(const blitz::Array<double, 2>& phasmoms_in) {
    phasmoms_ = phasmoms_in;
  }

  

  const blitz::Array<double, 2>& phasfunc_dn() const {
    return phasfunc_dn_;
  }

  void phasfunc_dn(const blitz::Array<double, 2>& phasfunc_dn_in) {
    phasfunc_dn_ = phasfunc_dn_in;
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

  

  const blitz::Array<double, 3>& l_phasmoms() const {
    return l_phasmoms_;
  }

  void l_phasmoms(const blitz::Array<double, 3>& l_phasmoms_in) {
    l_phasmoms_ = l_phasmoms_in;
  }

  

  const blitz::Array<double, 3>& l_phasfunc_dn() const {
    return l_phasfunc_dn_;
  }

  void l_phasfunc_dn(const blitz::Array<double, 3>& l_phasfunc_dn_in) {
    l_phasfunc_dn_ = l_phasfunc_dn_in;
  }

  

  const blitz::Array<double, 1>& mu1() const {
    return mu1_;
  }

  void mu1(const blitz::Array<double, 1>& mu1_in) {
    mu1_ = mu1_in;
  }

  

  const blitz::Array<double, 2>& legpoly_dn() const {
    return legpoly_dn_;
  }

  void legpoly_dn(const blitz::Array<double, 2>& legpoly_dn_in) {
    legpoly_dn_ = legpoly_dn_in;
  }

  

  const blitz::Array<double, 2>& losw_paths() const {
    return losw_paths_;
  }

  void losw_paths(const blitz::Array<double, 2>& losw_paths_in) {
    losw_paths_ = losw_paths_in;
  }

  

  const blitz::Array<double, 2>& losp_paths() const {
    return losp_paths_;
  }

  void losp_paths(const blitz::Array<double, 2>& losp_paths_in) {
    losp_paths_ = losp_paths_in;
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

  

  const blitz::Array<double, 3>& xfine_p() const {
    return xfine_p_;
  }

  void xfine_p(const blitz::Array<double, 3>& xfine_p_in) {
    xfine_p_ = xfine_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_p() const {
    return wfine_p_;
  }

  void wfine_p(const blitz::Array<double, 3>& wfine_p_in) {
    wfine_p_ = wfine_p_in;
  }

  

  const blitz::Array<double, 3>& sunpaths_p() const {
    return sunpaths_p_;
  }

  void sunpaths_p(const blitz::Array<double, 3>& sunpaths_p_in) {
    sunpaths_p_ = sunpaths_p_in;
  }

  

  const blitz::Array<int, 2>& ntraverse_p() const {
    return ntraverse_p_;
  }

  void ntraverse_p(const blitz::Array<int, 2>& ntraverse_p_in) {
    ntraverse_p_ = ntraverse_p_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine_p() const {
    return sunpathsfine_p_;
  }

  void sunpathsfine_p(const blitz::Array<double, 4>& sunpathsfine_p_in) {
    sunpathsfine_p_ = sunpathsfine_p_in;
  }

  

  const blitz::Array<int, 3>& ntraversefine_p() const {
    return ntraversefine_p_;
  }

  void ntraversefine_p(const blitz::Array<int, 3>& ntraversefine_p_in) {
    ntraversefine_p_ = ntraversefine_p_in;
  }

  

  const blitz::Array<double, 2>& intensity_dn() const {
    return intensity_dn_;
  }

  

  const blitz::Array<double, 4>& lp_jacobians_dn() const {
    return lp_jacobians_dn_;
  }

  

  const blitz::Array<double, 2>& lostrans_dn() const {
    return lostrans_dn_;
  }

  

  const blitz::Array<double, 3>& lp_lostrans_dn() const {
    return lp_lostrans_dn_;
  }

  

  const int& max_surfacewfs() const {
    return max_surfacewfs_;
  }

  

  const int& max_sleavewfs() const {
    return max_sleavewfs_;
  }

  

  const bool& do_surface_leaving() const {
    return do_surface_leaving_;
  }

  void do_surface_leaving(const bool& do_surface_leaving_in) {
    do_surface_leaving_ = do_surface_leaving_in;
  }

  

  const bool& do_water_leaving() const {
    return do_water_leaving_;
  }

  void do_water_leaving(const bool& do_water_leaving_in) {
    do_water_leaving_ = do_water_leaving_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_up() const {
    return do_sources_up_;
  }

  void do_sources_up(const blitz::Array<bool, 2>& do_sources_up_in) {
    do_sources_up_ = do_sources_up_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_up_p() const {
    return do_sources_up_p_;
  }

  void do_sources_up_p(const blitz::Array<bool, 2>& do_sources_up_p_in) {
    do_sources_up_p_ = do_sources_up_p_in;
  }

  

  const bool& do_surfacewfs() const {
    return do_surfacewfs_;
  }

  void do_surfacewfs(const bool& do_surfacewfs_in) {
    do_surfacewfs_ = do_surfacewfs_in;
  }

  

  const bool& do_sleavewfs() const {
    return do_sleavewfs_;
  }

  void do_sleavewfs(const bool& do_sleavewfs_in) {
    do_sleavewfs_ = do_sleavewfs_in;
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

  

  const int& n_surfacewfs() const {
    return n_surfacewfs_;
  }

  

  const blitz::Array<double, 2>& phasfunc_up() const {
    return phasfunc_up_;
  }

  void phasfunc_up(const blitz::Array<double, 2>& phasfunc_up_in) {
    phasfunc_up_ = phasfunc_up_in;
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

  

  const blitz::Array<double, 3>& l_phasfunc_up() const {
    return l_phasfunc_up_;
  }

  void l_phasfunc_up(const blitz::Array<double, 3>& l_phasfunc_up_in) {
    l_phasfunc_up_ = l_phasfunc_up_in;
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

  

  const blitz::Array<double, 1>& mu0() const {
    return mu0_;
  }

  void mu0(const blitz::Array<double, 1>& mu0_in) {
    mu0_ = mu0_in;
  }

  

  const blitz::Array<double, 2>& legpoly_up() const {
    return legpoly_up_;
  }

  void legpoly_up(const blitz::Array<double, 2>& legpoly_up_in) {
    legpoly_up_ = legpoly_up_in;
  }

  

  const blitz::Array<double, 2>& intensity_up() const {
    return intensity_up_;
  }

  

  const blitz::Array<double, 2>& intensity_db() const {
    return intensity_db_;
  }

  

  const blitz::Array<double, 4>& lp_jacobians_up() const {
    return lp_jacobians_up_;
  }

  

  const blitz::Array<double, 4>& lp_jacobians_db() const {
    return lp_jacobians_db_;
  }

  

  const blitz::Array<double, 3>& ls_jacobians_db() const {
    return ls_jacobians_db_;
  }

  

  const blitz::Array<double, 2>& cumtrans() const {
    return cumtrans_;
  }

  

  const blitz::Array<double, 2>& lostrans_up() const {
    return lostrans_up_;
  }

  

  const blitz::Array<double, 4>& lp_cumtrans() const {
    return lp_cumtrans_;
  }

  

  const blitz::Array<double, 3>& lp_lostrans_up() const {
    return lp_lostrans_up_;
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

  

  const blitz::Array<double, 3>& xfine_up() const {
    return xfine_up_;
  }

  void xfine_up(const blitz::Array<double, 3>& xfine_up_in) {
    xfine_up_ = xfine_up_in;
  }

  

  const blitz::Array<double, 3>& wfine_up() const {
    return wfine_up_;
  }

  void wfine_up(const blitz::Array<double, 3>& wfine_up_in) {
    wfine_up_ = wfine_up_in;
  }

  

  const blitz::Array<double, 3>& sunpaths_up() const {
    return sunpaths_up_;
  }

  void sunpaths_up(const blitz::Array<double, 3>& sunpaths_up_in) {
    sunpaths_up_ = sunpaths_up_in;
  }

  

  const blitz::Array<int, 2>& ntraverse_up() const {
    return ntraverse_up_;
  }

  void ntraverse_up(const blitz::Array<int, 2>& ntraverse_up_in) {
    ntraverse_up_ = ntraverse_up_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine_up() const {
    return sunpathsfine_up_;
  }

  void sunpathsfine_up(const blitz::Array<double, 4>& sunpathsfine_up_in) {
    sunpathsfine_up_ = sunpathsfine_up_in;
  }

  

  const blitz::Array<int, 3>& ntraversefine_up() const {
    return ntraversefine_up_;
  }

  void ntraversefine_up(const blitz::Array<int, 3>& ntraversefine_up_in) {
    ntraversefine_up_ = ntraversefine_up_in;
  }

  

  const blitz::Array<double, 3>& xfine_dn() const {
    return xfine_dn_;
  }

  void xfine_dn(const blitz::Array<double, 3>& xfine_dn_in) {
    xfine_dn_ = xfine_dn_in;
  }

  

  const blitz::Array<double, 3>& wfine_dn() const {
    return wfine_dn_;
  }

  void wfine_dn(const blitz::Array<double, 3>& wfine_dn_in) {
    wfine_dn_ = wfine_dn_in;
  }

  

  const blitz::Array<double, 3>& sunpaths_dn() const {
    return sunpaths_dn_;
  }

  void sunpaths_dn(const blitz::Array<double, 3>& sunpaths_dn_in) {
    sunpaths_dn_ = sunpaths_dn_in;
  }

  

  const blitz::Array<int, 2>& ntraverse_dn() const {
    return ntraverse_dn_;
  }

  void ntraverse_dn(const blitz::Array<int, 2>& ntraverse_dn_in) {
    ntraverse_dn_ = ntraverse_dn_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine_dn() const {
    return sunpathsfine_dn_;
  }

  void sunpathsfine_dn(const blitz::Array<double, 4>& sunpathsfine_dn_in) {
    sunpathsfine_dn_ = sunpathsfine_dn_in;
  }

  

  const blitz::Array<int, 3>& ntraversefine_dn() const {
    return ntraversefine_dn_;
  }

  void ntraversefine_dn(const blitz::Array<int, 3>& ntraversefine_dn_in) {
    ntraversefine_dn_ = ntraversefine_dn_in;
  }

  

  const blitz::Array<double, 3>& xfine_up_p() const {
    return xfine_up_p_;
  }

  void xfine_up_p(const blitz::Array<double, 3>& xfine_up_p_in) {
    xfine_up_p_ = xfine_up_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_up_p() const {
    return wfine_up_p_;
  }

  void wfine_up_p(const blitz::Array<double, 3>& wfine_up_p_in) {
    wfine_up_p_ = wfine_up_p_in;
  }

  

  const blitz::Array<double, 3>& sunpaths_up_p() const {
    return sunpaths_up_p_;
  }

  void sunpaths_up_p(const blitz::Array<double, 3>& sunpaths_up_p_in) {
    sunpaths_up_p_ = sunpaths_up_p_in;
  }

  

  const blitz::Array<int, 2>& ntraverse_up_p() const {
    return ntraverse_up_p_;
  }

  void ntraverse_up_p(const blitz::Array<int, 2>& ntraverse_up_p_in) {
    ntraverse_up_p_ = ntraverse_up_p_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine_up_p() const {
    return sunpathsfine_up_p_;
  }

  void sunpathsfine_up_p(const blitz::Array<double, 4>& sunpathsfine_up_p_in) {
    sunpathsfine_up_p_ = sunpathsfine_up_p_in;
  }

  

  const blitz::Array<int, 3>& ntraversefine_up_p() const {
    return ntraversefine_up_p_;
  }

  void ntraversefine_up_p(const blitz::Array<int, 3>& ntraversefine_up_p_in) {
    ntraversefine_up_p_ = ntraversefine_up_p_in;
  }

  

  const blitz::Array<double, 3>& xfine_dn_p() const {
    return xfine_dn_p_;
  }

  void xfine_dn_p(const blitz::Array<double, 3>& xfine_dn_p_in) {
    xfine_dn_p_ = xfine_dn_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_dn_p() const {
    return wfine_dn_p_;
  }

  void wfine_dn_p(const blitz::Array<double, 3>& wfine_dn_p_in) {
    wfine_dn_p_ = wfine_dn_p_in;
  }

  

  const blitz::Array<double, 3>& sunpaths_dn_p() const {
    return sunpaths_dn_p_;
  }

  void sunpaths_dn_p(const blitz::Array<double, 3>& sunpaths_dn_p_in) {
    sunpaths_dn_p_ = sunpaths_dn_p_in;
  }

  

  const blitz::Array<int, 2>& ntraverse_dn_p() const {
    return ntraverse_dn_p_;
  }

  void ntraverse_dn_p(const blitz::Array<int, 2>& ntraverse_dn_p_in) {
    ntraverse_dn_p_ = ntraverse_dn_p_in;
  }

  

  const blitz::Array<double, 4>& sunpathsfine_dn_p() const {
    return sunpathsfine_dn_p_;
  }

  void sunpathsfine_dn_p(const blitz::Array<double, 4>& sunpathsfine_dn_p_in) {
    sunpathsfine_dn_p_ = sunpathsfine_dn_p_in;
  }

  

  const blitz::Array<int, 3>& ntraversefine_dn_p() const {
    return ntraversefine_dn_p_;
  }

  void ntraversefine_dn_p(const blitz::Array<int, 3>& ntraversefine_dn_p_in) {
    ntraversefine_dn_p_ = ntraversefine_dn_p_in;
  }

  

  
  void ss_integral_ilps_dn() {
    
    
    fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_dn_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &maxmoments_input_, &max_user_levels_, &max_atmoswfs_, &do_deltam_scaling_, &do_phasfunc_, &do_partials_, &do_planpar_, &do_enhanced_ps_, &flux_, do_sources_dn_.dataFirst(), do_sources_dn_p_.dataFirst(), &do_profilewfs_, lvaryflags_.dataFirst(), lvarynums_.dataFirst(), lvarymoms_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &nmoments_input_, &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), phasmoms_.dataFirst(), phasfunc_dn_.dataFirst(), l_extinction_.dataFirst(), l_deltaus_.dataFirst(), l_omega_.dataFirst(), l_truncfac_.dataFirst(), l_phasmoms_.dataFirst(), l_phasfunc_dn_.dataFirst(), mu1_.dataFirst(), legpoly_dn_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), sunpaths_.dataFirst(), ntraverse_.dataFirst(), sunpathsfine_.dataFirst(), ntraversefine_.dataFirst(), xfine_p_.dataFirst(), wfine_p_.dataFirst(), sunpaths_p_.dataFirst(), ntraverse_p_.dataFirst(), sunpathsfine_p_.dataFirst(), ntraversefine_p_.dataFirst(), intensity_dn_.dataFirst(), lp_jacobians_dn_.dataFirst(), lostrans_dn_.dataFirst(), lp_lostrans_dn_.dataFirst());
    
  }
void ss_integral_ilps_up() {
    
    
    fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_up_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &maxmoments_input_, &max_user_levels_, &max_atmoswfs_, &max_surfacewfs_, &max_sleavewfs_, &do_deltam_scaling_, &do_phasfunc_, &do_surface_leaving_, &do_water_leaving_, &do_partials_, &do_planpar_, &do_enhanced_ps_, &flux_, do_sources_up_.dataFirst(), do_sources_up_p_.dataFirst(), &do_profilewfs_, &do_surfacewfs_, &do_sleavewfs_, &n_reflecwfs_, &n_sleavewfs_, &n_surfacewfs_, lvaryflags_.dataFirst(), lvarynums_.dataFirst(), lvarymoms_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &nmoments_input_, &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), phasmoms_.dataFirst(), phasfunc_up_.dataFirst(), reflec_.dataFirst(), slterm_.dataFirst(), l_extinction_.dataFirst(), l_deltaus_.dataFirst(), l_omega_.dataFirst(), l_truncfac_.dataFirst(), l_phasmoms_.dataFirst(), l_phasfunc_up_.dataFirst(), ls_reflec_.dataFirst(), lssl_slterm_.dataFirst(), mu0_.dataFirst(), mu1_.dataFirst(), legpoly_up_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), sunpaths_.dataFirst(), ntraverse_.dataFirst(), sunpathsfine_.dataFirst(), ntraversefine_.dataFirst(), xfine_p_.dataFirst(), wfine_p_.dataFirst(), sunpaths_p_.dataFirst(), ntraverse_p_.dataFirst(), sunpathsfine_p_.dataFirst(), ntraversefine_p_.dataFirst(), intensity_up_.dataFirst(), intensity_db_.dataFirst(), lp_jacobians_up_.dataFirst(), lp_jacobians_db_.dataFirst(), ls_jacobians_db_.dataFirst(), cumtrans_.dataFirst(), lostrans_up_.dataFirst(), lp_cumtrans_.dataFirst(), lp_lostrans_up_.dataFirst());
    
  }
void ss_integral_ilps_updn() {
    
    
    fo_scalarss_rtcalcs_ilps_m_ss_integral_ilps_updn_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &maxmoments_input_, &max_user_levels_, &max_atmoswfs_, &max_surfacewfs_, &max_sleavewfs_, &do_upwelling_, &do_dnwelling_, &do_deltam_scaling_, &do_phasfunc_, &do_surface_leaving_, &do_water_leaving_, &do_partials_, &do_planpar_, &do_enhanced_ps_, do_sources_up_.dataFirst(), do_sources_up_p_.dataFirst(), do_sources_dn_.dataFirst(), do_sources_dn_p_.dataFirst(), &do_profilewfs_, &do_surfacewfs_, &do_sleavewfs_, &n_reflecwfs_, &n_sleavewfs_, &n_surfacewfs_, lvaryflags_.dataFirst(), lvarynums_.dataFirst(), lvarymoms_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &nmoments_input_, &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), &flux_, extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), phasmoms_.dataFirst(), phasfunc_up_.dataFirst(), phasfunc_dn_.dataFirst(), reflec_.dataFirst(), slterm_.dataFirst(), ls_reflec_.dataFirst(), lssl_slterm_.dataFirst(), l_extinction_.dataFirst(), l_deltaus_.dataFirst(), l_omega_.dataFirst(), l_truncfac_.dataFirst(), l_phasmoms_.dataFirst(), l_phasfunc_up_.dataFirst(), l_phasfunc_dn_.dataFirst(), mu0_.dataFirst(), mu1_.dataFirst(), legpoly_up_.dataFirst(), legpoly_dn_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_up_.dataFirst(), wfine_up_.dataFirst(), sunpaths_up_.dataFirst(), ntraverse_up_.dataFirst(), sunpathsfine_up_.dataFirst(), ntraversefine_up_.dataFirst(), xfine_dn_.dataFirst(), wfine_dn_.dataFirst(), sunpaths_dn_.dataFirst(), ntraverse_dn_.dataFirst(), sunpathsfine_dn_.dataFirst(), ntraversefine_dn_.dataFirst(), xfine_up_p_.dataFirst(), wfine_up_p_.dataFirst(), sunpaths_up_p_.dataFirst(), ntraverse_up_p_.dataFirst(), sunpathsfine_up_p_.dataFirst(), ntraversefine_up_p_.dataFirst(), xfine_dn_p_.dataFirst(), wfine_dn_p_.dataFirst(), sunpaths_dn_p_.dataFirst(), ntraverse_dn_p_.dataFirst(), sunpathsfine_dn_p_.dataFirst(), ntraversefine_dn_p_.dataFirst(), intensity_up_.dataFirst(), intensity_db_.dataFirst(), lp_jacobians_up_.dataFirst(), lp_jacobians_db_.dataFirst(), ls_jacobians_db_.dataFirst(), lostrans_up_.dataFirst(), lp_lostrans_up_.dataFirst(), intensity_dn_.dataFirst(), lp_jacobians_dn_.dataFirst(), cumtrans_.dataFirst(), lp_cumtrans_.dataFirst(), lostrans_dn_.dataFirst(), lp_lostrans_dn_.dataFirst());
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "Fo_Scalarss_Rtcalcs_Ilps:" << std::endl
      << "          maxgeoms: " << maxgeoms()  << std::endl
      << "         maxlayers: " << maxlayers()  << std::endl
      << "       maxpartials: " << maxpartials()  << std::endl
      << "           maxfine: " << maxfine()  << std::endl
      << "  maxmoments_input: " << maxmoments_input()  << std::endl
      << "   max_user_levels: " << max_user_levels()  << std::endl
      << "      max_atmoswfs: " << max_atmoswfs()  << std::endl
      << " do_deltam_scaling: " << do_deltam_scaling()  << std::endl
      << "       do_phasfunc: " << do_phasfunc()  << std::endl
      << "       do_partials: " << do_partials()  << std::endl
      << "        do_planpar: " << do_planpar()  << std::endl
      << "    do_enhanced_ps: " << do_enhanced_ps()  << std::endl
      << "              flux: " << flux()  << std::endl
      << "     do_sources_dn: " << std::endl << do_sources_dn()  << std::endl
      << "   do_sources_dn_p: " << std::endl << do_sources_dn_p()  << std::endl
      << "     do_profilewfs: " << do_profilewfs()  << std::endl
      << "        lvaryflags: " << std::endl << lvaryflags()  << std::endl
      << "         lvarynums: " << std::endl << lvarynums()  << std::endl
      << "         lvarymoms: " << std::endl << lvarymoms()  << std::endl
      << "            ngeoms: " << ngeoms()  << std::endl
      << "           nlayers: " << nlayers()  << std::endl
      << "         nfinedivs: " << std::endl << nfinedivs()  << std::endl
      << "    nmoments_input: " << nmoments_input()  << std::endl
      << "     n_user_levels: " << n_user_levels()  << std::endl
      << "       user_levels: " << std::endl << user_levels()  << std::endl
      << "         npartials: " << npartials()  << std::endl
      << "       nfinedivs_p: " << std::endl << nfinedivs_p()  << std::endl
      << "  partial_outindex: " << std::endl << partial_outindex()  << std::endl
      << "   partial_outflag: " << std::endl << partial_outflag()  << std::endl
      << "  partial_layeridx: " << std::endl << partial_layeridx()  << std::endl
      << "        extinction: " << std::endl << extinction()  << std::endl
      << "           deltaus: " << std::endl << deltaus()  << std::endl
      << "             omega: " << std::endl << omega()  << std::endl
      << "          truncfac: " << std::endl << truncfac()  << std::endl
      << "          phasmoms: " << std::endl << phasmoms()  << std::endl
      << "       phasfunc_dn: " << std::endl << phasfunc_dn()  << std::endl
      << "      l_extinction: " << std::endl << l_extinction()  << std::endl
      << "         l_deltaus: " << std::endl << l_deltaus()  << std::endl
      << "           l_omega: " << std::endl << l_omega()  << std::endl
      << "        l_truncfac: " << std::endl << l_truncfac()  << std::endl
      << "        l_phasmoms: " << std::endl << l_phasmoms()  << std::endl
      << "     l_phasfunc_dn: " << std::endl << l_phasfunc_dn()  << std::endl
      << "               mu1: " << std::endl << mu1()  << std::endl
      << "        legpoly_dn: " << std::endl << legpoly_dn()  << std::endl
      << "        losw_paths: " << std::endl << losw_paths()  << std::endl
      << "        losp_paths: " << std::endl << losp_paths()  << std::endl
      << "             xfine: " << std::endl << xfine()  << std::endl
      << "             wfine: " << std::endl << wfine()  << std::endl
      << "          sunpaths: " << std::endl << sunpaths()  << std::endl
      << "         ntraverse: " << std::endl << ntraverse()  << std::endl
      << "      sunpathsfine: " << std::endl << sunpathsfine()  << std::endl
      << "     ntraversefine: " << std::endl << ntraversefine()  << std::endl
      << "           xfine_p: " << std::endl << xfine_p()  << std::endl
      << "           wfine_p: " << std::endl << wfine_p()  << std::endl
      << "        sunpaths_p: " << std::endl << sunpaths_p()  << std::endl
      << "       ntraverse_p: " << std::endl << ntraverse_p()  << std::endl
      << "    sunpathsfine_p: " << std::endl << sunpathsfine_p()  << std::endl
      << "   ntraversefine_p: " << std::endl << ntraversefine_p()  << std::endl
      << "      intensity_dn: " << std::endl << intensity_dn()  << std::endl
      << "   lp_jacobians_dn: " << std::endl << lp_jacobians_dn()  << std::endl
      << "       lostrans_dn: " << std::endl << lostrans_dn()  << std::endl
      << "    lp_lostrans_dn: " << std::endl << lp_lostrans_dn()  << std::endl
      << "    max_surfacewfs: " << max_surfacewfs()  << std::endl
      << "     max_sleavewfs: " << max_sleavewfs()  << std::endl
      << "do_surface_leaving: " << do_surface_leaving()  << std::endl
      << "  do_water_leaving: " << do_water_leaving()  << std::endl
      << "     do_sources_up: " << std::endl << do_sources_up()  << std::endl
      << "   do_sources_up_p: " << std::endl << do_sources_up_p()  << std::endl
      << "     do_surfacewfs: " << do_surfacewfs()  << std::endl
      << "      do_sleavewfs: " << do_sleavewfs()  << std::endl
      << "       n_reflecwfs: " << n_reflecwfs()  << std::endl
      << "       n_sleavewfs: " << n_sleavewfs()  << std::endl
      << "      n_surfacewfs: " << n_surfacewfs()  << std::endl
      << "       phasfunc_up: " << std::endl << phasfunc_up()  << std::endl
      << "            reflec: " << std::endl << reflec()  << std::endl
      << "            slterm: " << std::endl << slterm()  << std::endl
      << "     l_phasfunc_up: " << std::endl << l_phasfunc_up()  << std::endl
      << "         ls_reflec: " << std::endl << ls_reflec()  << std::endl
      << "       lssl_slterm: " << std::endl << lssl_slterm()  << std::endl
      << "               mu0: " << std::endl << mu0()  << std::endl
      << "        legpoly_up: " << std::endl << legpoly_up()  << std::endl
      << "      intensity_up: " << std::endl << intensity_up()  << std::endl
      << "      intensity_db: " << std::endl << intensity_db()  << std::endl
      << "   lp_jacobians_up: " << std::endl << lp_jacobians_up()  << std::endl
      << "   lp_jacobians_db: " << std::endl << lp_jacobians_db()  << std::endl
      << "   ls_jacobians_db: " << std::endl << ls_jacobians_db()  << std::endl
      << "          cumtrans: " << std::endl << cumtrans()  << std::endl
      << "       lostrans_up: " << std::endl << lostrans_up()  << std::endl
      << "       lp_cumtrans: " << std::endl << lp_cumtrans()  << std::endl
      << "    lp_lostrans_up: " << std::endl << lp_lostrans_up()  << std::endl
      << "      do_upwelling: " << do_upwelling()  << std::endl
      << "      do_dnwelling: " << do_dnwelling()  << std::endl
      << "          xfine_up: " << std::endl << xfine_up()  << std::endl
      << "          wfine_up: " << std::endl << wfine_up()  << std::endl
      << "       sunpaths_up: " << std::endl << sunpaths_up()  << std::endl
      << "      ntraverse_up: " << std::endl << ntraverse_up()  << std::endl
      << "   sunpathsfine_up: " << std::endl << sunpathsfine_up()  << std::endl
      << "  ntraversefine_up: " << std::endl << ntraversefine_up()  << std::endl
      << "          xfine_dn: " << std::endl << xfine_dn()  << std::endl
      << "          wfine_dn: " << std::endl << wfine_dn()  << std::endl
      << "       sunpaths_dn: " << std::endl << sunpaths_dn()  << std::endl
      << "      ntraverse_dn: " << std::endl << ntraverse_dn()  << std::endl
      << "   sunpathsfine_dn: " << std::endl << sunpathsfine_dn()  << std::endl
      << "  ntraversefine_dn: " << std::endl << ntraversefine_dn()  << std::endl
      << "        xfine_up_p: " << std::endl << xfine_up_p()  << std::endl
      << "        wfine_up_p: " << std::endl << wfine_up_p()  << std::endl
      << "     sunpaths_up_p: " << std::endl << sunpaths_up_p()  << std::endl
      << "    ntraverse_up_p: " << std::endl << ntraverse_up_p()  << std::endl
      << " sunpathsfine_up_p: " << std::endl << sunpathsfine_up_p()  << std::endl
      << "ntraversefine_up_p: " << std::endl << ntraversefine_up_p()  << std::endl
      << "        xfine_dn_p: " << std::endl << xfine_dn_p()  << std::endl
      << "        wfine_dn_p: " << std::endl << wfine_dn_p()  << std::endl
      << "     sunpaths_dn_p: " << std::endl << sunpaths_dn_p()  << std::endl
      << "    ntraverse_dn_p: " << std::endl << ntraverse_dn_p()  << std::endl
      << " sunpathsfine_dn_p: " << std::endl << sunpathsfine_dn_p()  << std::endl
      << "ntraversefine_dn_p: " << std::endl << ntraversefine_dn_p()  << std::endl;

  }

private:
  int maxgeoms_;
  int maxlayers_;
  int maxpartials_;
  int maxfine_;
  int maxmoments_input_;
  int max_user_levels_;
  int max_atmoswfs_;
  bool do_deltam_scaling_;
  bool do_phasfunc_;
  bool do_partials_;
  bool do_planpar_;
  bool do_enhanced_ps_;
  double flux_;
  blitz::Array<bool, 2> do_sources_dn_;
  blitz::Array<bool, 2> do_sources_dn_p_;
  bool do_profilewfs_;
  blitz::Array<bool, 1> lvaryflags_;
  blitz::Array<int, 1> lvarynums_;
  blitz::Array<bool, 2> lvarymoms_;
  int ngeoms_;
  int nlayers_;
  blitz::Array<int, 2> nfinedivs_;
  int nmoments_input_;
  int n_user_levels_;
  blitz::Array<int, 1> user_levels_;
  int npartials_;
  blitz::Array<int, 2> nfinedivs_p_;
  blitz::Array<int, 1> partial_outindex_;
  blitz::Array<bool, 1> partial_outflag_;
  blitz::Array<int, 1> partial_layeridx_;
  blitz::Array<double, 1> extinction_;
  blitz::Array<double, 1> deltaus_;
  blitz::Array<double, 1> omega_;
  blitz::Array<double, 1> truncfac_;
  blitz::Array<double, 2> phasmoms_;
  blitz::Array<double, 2> phasfunc_dn_;
  blitz::Array<double, 2> l_extinction_;
  blitz::Array<double, 2> l_deltaus_;
  blitz::Array<double, 2> l_omega_;
  blitz::Array<double, 2> l_truncfac_;
  blitz::Array<double, 3> l_phasmoms_;
  blitz::Array<double, 3> l_phasfunc_dn_;
  blitz::Array<double, 1> mu1_;
  blitz::Array<double, 2> legpoly_dn_;
  blitz::Array<double, 2> losw_paths_;
  blitz::Array<double, 2> losp_paths_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> sunpaths_;
  blitz::Array<int, 2> ntraverse_;
  blitz::Array<double, 4> sunpathsfine_;
  blitz::Array<int, 3> ntraversefine_;
  blitz::Array<double, 3> xfine_p_;
  blitz::Array<double, 3> wfine_p_;
  blitz::Array<double, 3> sunpaths_p_;
  blitz::Array<int, 2> ntraverse_p_;
  blitz::Array<double, 4> sunpathsfine_p_;
  blitz::Array<int, 3> ntraversefine_p_;
  blitz::Array<double, 2> intensity_dn_;
  blitz::Array<double, 4> lp_jacobians_dn_;
  blitz::Array<double, 2> lostrans_dn_;
  blitz::Array<double, 3> lp_lostrans_dn_;
  int max_surfacewfs_;
  int max_sleavewfs_;
  bool do_surface_leaving_;
  bool do_water_leaving_;
  blitz::Array<bool, 2> do_sources_up_;
  blitz::Array<bool, 2> do_sources_up_p_;
  bool do_surfacewfs_;
  bool do_sleavewfs_;
  int n_reflecwfs_;
  int n_sleavewfs_;
  int n_surfacewfs_;
  blitz::Array<double, 2> phasfunc_up_;
  blitz::Array<double, 1> reflec_;
  blitz::Array<double, 1> slterm_;
  blitz::Array<double, 3> l_phasfunc_up_;
  blitz::Array<double, 2> ls_reflec_;
  blitz::Array<double, 2> lssl_slterm_;
  blitz::Array<double, 1> mu0_;
  blitz::Array<double, 2> legpoly_up_;
  blitz::Array<double, 2> intensity_up_;
  blitz::Array<double, 2> intensity_db_;
  blitz::Array<double, 4> lp_jacobians_up_;
  blitz::Array<double, 4> lp_jacobians_db_;
  blitz::Array<double, 3> ls_jacobians_db_;
  blitz::Array<double, 2> cumtrans_;
  blitz::Array<double, 2> lostrans_up_;
  blitz::Array<double, 4> lp_cumtrans_;
  blitz::Array<double, 3> lp_lostrans_up_;
  bool do_upwelling_;
  bool do_dnwelling_;
  blitz::Array<double, 3> xfine_up_;
  blitz::Array<double, 3> wfine_up_;
  blitz::Array<double, 3> sunpaths_up_;
  blitz::Array<int, 2> ntraverse_up_;
  blitz::Array<double, 4> sunpathsfine_up_;
  blitz::Array<int, 3> ntraversefine_up_;
  blitz::Array<double, 3> xfine_dn_;
  blitz::Array<double, 3> wfine_dn_;
  blitz::Array<double, 3> sunpaths_dn_;
  blitz::Array<int, 2> ntraverse_dn_;
  blitz::Array<double, 4> sunpathsfine_dn_;
  blitz::Array<int, 3> ntraversefine_dn_;
  blitz::Array<double, 3> xfine_up_p_;
  blitz::Array<double, 3> wfine_up_p_;
  blitz::Array<double, 3> sunpaths_up_p_;
  blitz::Array<int, 2> ntraverse_up_p_;
  blitz::Array<double, 4> sunpathsfine_up_p_;
  blitz::Array<int, 3> ntraversefine_up_p_;
  blitz::Array<double, 3> xfine_dn_p_;
  blitz::Array<double, 3> wfine_dn_p_;
  blitz::Array<double, 3> sunpaths_dn_p_;
  blitz::Array<int, 2> ntraverse_dn_p_;
  blitz::Array<double, 4> sunpathsfine_dn_p_;
  blitz::Array<int, 3> ntraversefine_dn_p_;

  // Serialization support
  Fo_Scalarss_Rtcalcs_Ilps() {}
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "fo_scalarss_spherfuncs_m" in file: "FO_ScalarSS_Spherfuncs.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_scalarss_spherfuncs_m_fo_scalarss_spherfuncs_wrap(const int* maxmoments_in, const int* maxgeoms_in, const int* nmoments_in, const int* ngeoms_in, const bool* starter_in, const bool* do_spherfunc_in, const double* cosscat_in, const double* df1_in, const double* df2_in, const double* ss_pleg_in);
}

class Fo_Scalarss_Spherfuncs : public Printable<Fo_Scalarss_Spherfuncs> {

public:
  Fo_Scalarss_Spherfuncs(const int& maxmoments_in, const int& maxgeoms_in, const int& nmoments_in, const int& ngeoms_in) : maxmoments_(maxmoments_in), maxgeoms_(maxgeoms_in), nmoments_(nmoments_in), ngeoms_(ngeoms_in) 
  { 
    starter_ = false;
    do_spherfunc_ = false;
    cosscat_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    cosscat_ = 0;
    df1_.reference( blitz::Array<double, 1>(maxmoments_, blitz::ColumnMajorArray<1>()) );
    df1_ = 0;
    df2_.reference( blitz::Array<double, 1>(maxmoments_, blitz::ColumnMajorArray<1>()) );
    df2_ = 0;
    ss_pleg_.reference( blitz::Array<double, 2>(maxmoments_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    ss_pleg_ = 0;
    // Initialize type pointers
    
  }

  virtual ~Fo_Scalarss_Spherfuncs() = default;

  const int& maxmoments() const {
    return maxmoments_;
  }

  

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& nmoments() const {
    return nmoments_;
  }

  

  const int& ngeoms() const {
    return ngeoms_;
  }

  

  const bool& starter() const {
    return starter_;
  }

  void starter(const bool& starter_in) {
    starter_ = starter_in;
  }

  

  const bool& do_spherfunc() const {
    return do_spherfunc_;
  }

  void do_spherfunc(const bool& do_spherfunc_in) {
    do_spherfunc_ = do_spherfunc_in;
  }

  

  const blitz::Array<double, 1>& cosscat() const {
    return cosscat_;
  }

  void cosscat(const blitz::Array<double, 1>& cosscat_in) {
    cosscat_ = cosscat_in;
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

  

  const blitz::Array<double, 2>& ss_pleg() const {
    return ss_pleg_;
  }

  

  
  void run() {
    
    
    fo_scalarss_spherfuncs_m_fo_scalarss_spherfuncs_wrap(&maxmoments_, &maxgeoms_, &nmoments_, &ngeoms_, &starter_, &do_spherfunc_, cosscat_.dataFirst(), df1_.dataFirst(), df2_.dataFirst(), ss_pleg_.dataFirst());
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "Fo_Scalarss_Spherfuncs:" << std::endl
      << "  maxmoments: " << maxmoments()  << std::endl
      << "    maxgeoms: " << maxgeoms()  << std::endl
      << "    nmoments: " << nmoments()  << std::endl
      << "      ngeoms: " << ngeoms()  << std::endl
      << "     starter: " << starter()  << std::endl
      << "do_spherfunc: " << do_spherfunc()  << std::endl
      << "     cosscat: " << std::endl << cosscat()  << std::endl
      << "         df1: " << std::endl << df1()  << std::endl
      << "         df2: " << std::endl << df2()  << std::endl
      << "     ss_pleg: " << std::endl << ss_pleg()  << std::endl;

  }

private:
  int maxmoments_;
  int maxgeoms_;
  int nmoments_;
  int ngeoms_;
  bool starter_;
  bool do_spherfunc_;
  blitz::Array<double, 1> cosscat_;
  blitz::Array<double, 1> df1_;
  blitz::Array<double, 1> df2_;
  blitz::Array<double, 2> ss_pleg_;

  // Serialization support
  Fo_Scalarss_Spherfuncs() {}
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "fo_thermal_rtcalcs_i_m" in file: "FO_Thermal_RTCalcs_I.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_thermal_rtcalcs_i_m_dte_integral_i_dn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* max_user_levels_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const bool* do_sources_dn_in, const bool* do_sources_dn_p_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* bb_input_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* mu1_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_in, const double* wfine_in, const double* hfine_in, const double* xfine_p_in, const double* wfine_p_in, const double* hfine_p_in, const double* intensity_dta_dn_in, const double* cumsource_dn_in, const double* tcom1_in);
  void fo_thermal_rtcalcs_i_m_dte_integral_i_up_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* max_user_levels_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const bool* do_sources_up_in, const bool* do_sources_up_p_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* bb_input_in, const double* surfbb_in, const double* user_emissivity_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* mu1_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_in, const double* wfine_in, const double* hfine_in, const double* xfine_p_in, const double* wfine_p_in, const double* hfine_p_in, const double* intensity_dta_up_in, const double* intensity_dts_in, const double* cumsource_up_in, const double* tcom1_in, const double* lostrans_up_in, const double* lostrans_up_p_in);
  void fo_thermal_rtcalcs_i_m_dte_integral_i_updn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* max_user_levels_in, const bool* do_upwelling_in, const bool* do_dnwelling_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const bool* do_sources_up_in, const bool* do_sources_up_p_in, const bool* do_sources_dn_in, const bool* do_sources_dn_p_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* bb_input_in, const double* surfbb_in, const double* user_emissivity_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* mu1_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_up_in, const double* wfine_up_in, const double* hfine_up_in, const double* xfine_dn_in, const double* wfine_dn_in, const double* hfine_dn_in, const double* xfine_up_p_in, const double* wfine_up_p_in, const double* hfine_up_p_in, const double* xfine_dn_p_in, const double* wfine_dn_p_in, const double* hfine_dn_p_in, const double* intensity_dta_up_in, const double* intensity_dts_in, const double* intensity_dta_dn_in, const double* cumsource_up_in, const double* cumsource_dn_in, const double* tcom1_in, const double* lostrans_up_in, const double* lostrans_up_p_in);
}

class Fo_Thermal_Rtcalcs_I : public Printable<Fo_Thermal_Rtcalcs_I> {

public:
  Fo_Thermal_Rtcalcs_I(const int& maxgeoms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& max_user_levels_in, const int& ngeoms_in, const int& nlayers_in, const int& n_user_levels_in, const int& npartials_in) : maxgeoms_(maxgeoms_in), maxlayers_(maxlayers_in), maxpartials_(maxpartials_in), maxfine_(maxfine_in), max_user_levels_(max_user_levels_in), ngeoms_(ngeoms_in), nlayers_(nlayers_in), n_user_levels_(n_user_levels_in), npartials_(npartials_in) 
  { 
    do_thermset_ = false;
    do_deltam_scaling_ = false;
    do_partials_ = false;
    do_planpar_ = false;
    do_enhanced_ps_ = false;
    do_sources_dn_.reference( blitz::Array<bool, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_dn_ = false;
    do_sources_dn_p_.reference( blitz::Array<bool, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_dn_p_ = false;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    user_levels_.reference( blitz::Array<int, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    user_levels_ = 0;
    nfinedivs_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_p_ = 0;
    partial_outindex_.reference( blitz::Array<int, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    partial_outindex_ = 0;
    partial_outflag_.reference( blitz::Array<bool, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    partial_outflag_ = false;
    partial_layeridx_.reference( blitz::Array<int, 1>(maxpartials_, blitz::ColumnMajorArray<1>()) );
    partial_layeridx_ = 0;
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
    losw_paths_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losw_paths_ = 0;
    losp_paths_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losp_paths_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    hfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_ = 0;
    xfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_p_ = 0;
    wfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_p_ = 0;
    hfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_p_ = 0;
    intensity_dta_dn_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dta_dn_ = 0;
    cumsource_dn_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cumsource_dn_ = 0;
    tcom1_.reference( blitz::Array<double, 2>(maxlayers_, 2, blitz::ColumnMajorArray<2>()) );
    tcom1_ = 0;
    do_sources_up_.reference( blitz::Array<bool, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_up_ = false;
    do_sources_up_p_.reference( blitz::Array<bool, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_up_p_ = false;
    surfbb_ = 0;
    user_emissivity_.reference( blitz::Array<double, 1>(maxgeoms_, blitz::ColumnMajorArray<1>()) );
    user_emissivity_ = 0;
    intensity_dta_up_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dta_up_ = 0;
    intensity_dts_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dts_ = 0;
    cumsource_up_.reference( blitz::Array<double, 2>(maxlayers_-0+1, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    cumsource_up_ = 0;
    lostrans_up_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    lostrans_up_ = 0;
    lostrans_up_p_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    lostrans_up_p_ = 0;
    do_upwelling_ = false;
    do_dnwelling_ = false;
    xfine_up_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_up_ = 0;
    wfine_up_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_up_ = 0;
    hfine_up_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_up_ = 0;
    xfine_dn_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_dn_ = 0;
    wfine_dn_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_dn_ = 0;
    hfine_dn_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_dn_ = 0;
    xfine_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_up_p_ = 0;
    wfine_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_up_p_ = 0;
    hfine_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_up_p_ = 0;
    xfine_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_dn_p_ = 0;
    wfine_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_dn_p_ = 0;
    hfine_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_dn_p_ = 0;
    // Initialize type pointers
    
  }

  virtual ~Fo_Thermal_Rtcalcs_I() = default;

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxpartials() const {
    return maxpartials_;
  }

  

  const int& maxfine() const {
    return maxfine_;
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

  

  const bool& do_partials() const {
    return do_partials_;
  }

  void do_partials(const bool& do_partials_in) {
    do_partials_ = do_partials_in;
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

  

  const blitz::Array<bool, 2>& do_sources_dn() const {
    return do_sources_dn_;
  }

  void do_sources_dn(const blitz::Array<bool, 2>& do_sources_dn_in) {
    do_sources_dn_ = do_sources_dn_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_dn_p() const {
    return do_sources_dn_p_;
  }

  void do_sources_dn_p(const blitz::Array<bool, 2>& do_sources_dn_p_in) {
    do_sources_dn_p_ = do_sources_dn_p_in;
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

  

  const int& npartials() const {
    return npartials_;
  }

  

  const blitz::Array<int, 2>& nfinedivs_p() const {
    return nfinedivs_p_;
  }

  void nfinedivs_p(const blitz::Array<int, 2>& nfinedivs_p_in) {
    nfinedivs_p_ = nfinedivs_p_in;
  }

  

  const blitz::Array<int, 1>& partial_outindex() const {
    return partial_outindex_;
  }

  void partial_outindex(const blitz::Array<int, 1>& partial_outindex_in) {
    partial_outindex_ = partial_outindex_in;
  }

  

  const blitz::Array<bool, 1>& partial_outflag() const {
    return partial_outflag_;
  }

  void partial_outflag(const blitz::Array<bool, 1>& partial_outflag_in) {
    partial_outflag_ = partial_outflag_in;
  }

  

  const blitz::Array<int, 1>& partial_layeridx() const {
    return partial_layeridx_;
  }

  void partial_layeridx(const blitz::Array<int, 1>& partial_layeridx_in) {
    partial_layeridx_ = partial_layeridx_in;
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

  

  const blitz::Array<double, 2>& losw_paths() const {
    return losw_paths_;
  }

  void losw_paths(const blitz::Array<double, 2>& losw_paths_in) {
    losw_paths_ = losw_paths_in;
  }

  

  const blitz::Array<double, 2>& losp_paths() const {
    return losp_paths_;
  }

  void losp_paths(const blitz::Array<double, 2>& losp_paths_in) {
    losp_paths_ = losp_paths_in;
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

  

  const blitz::Array<double, 3>& hfine() const {
    return hfine_;
  }

  void hfine(const blitz::Array<double, 3>& hfine_in) {
    hfine_ = hfine_in;
  }

  

  const blitz::Array<double, 3>& xfine_p() const {
    return xfine_p_;
  }

  void xfine_p(const blitz::Array<double, 3>& xfine_p_in) {
    xfine_p_ = xfine_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_p() const {
    return wfine_p_;
  }

  void wfine_p(const blitz::Array<double, 3>& wfine_p_in) {
    wfine_p_ = wfine_p_in;
  }

  

  const blitz::Array<double, 3>& hfine_p() const {
    return hfine_p_;
  }

  void hfine_p(const blitz::Array<double, 3>& hfine_p_in) {
    hfine_p_ = hfine_p_in;
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

  

  const blitz::Array<bool, 2>& do_sources_up() const {
    return do_sources_up_;
  }

  void do_sources_up(const blitz::Array<bool, 2>& do_sources_up_in) {
    do_sources_up_ = do_sources_up_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_up_p() const {
    return do_sources_up_p_;
  }

  void do_sources_up_p(const blitz::Array<bool, 2>& do_sources_up_p_in) {
    do_sources_up_p_ = do_sources_up_p_in;
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

  

  const blitz::Array<double, 2>& lostrans_up() const {
    return lostrans_up_;
  }

  

  const blitz::Array<double, 2>& lostrans_up_p() const {
    return lostrans_up_p_;
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

  

  const blitz::Array<double, 3>& xfine_up() const {
    return xfine_up_;
  }

  void xfine_up(const blitz::Array<double, 3>& xfine_up_in) {
    xfine_up_ = xfine_up_in;
  }

  

  const blitz::Array<double, 3>& wfine_up() const {
    return wfine_up_;
  }

  void wfine_up(const blitz::Array<double, 3>& wfine_up_in) {
    wfine_up_ = wfine_up_in;
  }

  

  const blitz::Array<double, 3>& hfine_up() const {
    return hfine_up_;
  }

  void hfine_up(const blitz::Array<double, 3>& hfine_up_in) {
    hfine_up_ = hfine_up_in;
  }

  

  const blitz::Array<double, 3>& xfine_dn() const {
    return xfine_dn_;
  }

  void xfine_dn(const blitz::Array<double, 3>& xfine_dn_in) {
    xfine_dn_ = xfine_dn_in;
  }

  

  const blitz::Array<double, 3>& wfine_dn() const {
    return wfine_dn_;
  }

  void wfine_dn(const blitz::Array<double, 3>& wfine_dn_in) {
    wfine_dn_ = wfine_dn_in;
  }

  

  const blitz::Array<double, 3>& hfine_dn() const {
    return hfine_dn_;
  }

  void hfine_dn(const blitz::Array<double, 3>& hfine_dn_in) {
    hfine_dn_ = hfine_dn_in;
  }

  

  const blitz::Array<double, 3>& xfine_up_p() const {
    return xfine_up_p_;
  }

  void xfine_up_p(const blitz::Array<double, 3>& xfine_up_p_in) {
    xfine_up_p_ = xfine_up_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_up_p() const {
    return wfine_up_p_;
  }

  void wfine_up_p(const blitz::Array<double, 3>& wfine_up_p_in) {
    wfine_up_p_ = wfine_up_p_in;
  }

  

  const blitz::Array<double, 3>& hfine_up_p() const {
    return hfine_up_p_;
  }

  void hfine_up_p(const blitz::Array<double, 3>& hfine_up_p_in) {
    hfine_up_p_ = hfine_up_p_in;
  }

  

  const blitz::Array<double, 3>& xfine_dn_p() const {
    return xfine_dn_p_;
  }

  void xfine_dn_p(const blitz::Array<double, 3>& xfine_dn_p_in) {
    xfine_dn_p_ = xfine_dn_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_dn_p() const {
    return wfine_dn_p_;
  }

  void wfine_dn_p(const blitz::Array<double, 3>& wfine_dn_p_in) {
    wfine_dn_p_ = wfine_dn_p_in;
  }

  

  const blitz::Array<double, 3>& hfine_dn_p() const {
    return hfine_dn_p_;
  }

  void hfine_dn_p(const blitz::Array<double, 3>& hfine_dn_p_in) {
    hfine_dn_p_ = hfine_dn_p_in;
  }

  

  
  void dte_integral_i_dn() {
    
    
    fo_thermal_rtcalcs_i_m_dte_integral_i_dn_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &max_user_levels_, &do_thermset_, &do_deltam_scaling_, &do_partials_, &do_planpar_, &do_enhanced_ps_, do_sources_dn_.dataFirst(), do_sources_dn_p_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), bb_input_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), mu1_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), hfine_.dataFirst(), xfine_p_.dataFirst(), wfine_p_.dataFirst(), hfine_p_.dataFirst(), intensity_dta_dn_.dataFirst(), cumsource_dn_.dataFirst(), tcom1_.dataFirst());
    
  }
void dte_integral_i_up() {
    
    
    fo_thermal_rtcalcs_i_m_dte_integral_i_up_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &max_user_levels_, &do_thermset_, &do_deltam_scaling_, &do_partials_, &do_planpar_, &do_enhanced_ps_, do_sources_up_.dataFirst(), do_sources_up_p_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), bb_input_.dataFirst(), &surfbb_, user_emissivity_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), mu1_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), hfine_.dataFirst(), xfine_p_.dataFirst(), wfine_p_.dataFirst(), hfine_p_.dataFirst(), intensity_dta_up_.dataFirst(), intensity_dts_.dataFirst(), cumsource_up_.dataFirst(), tcom1_.dataFirst(), lostrans_up_.dataFirst(), lostrans_up_p_.dataFirst());
    
  }
void dte_integral_i_updn() {
    
    
    fo_thermal_rtcalcs_i_m_dte_integral_i_updn_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &max_user_levels_, &do_upwelling_, &do_dnwelling_, &do_thermset_, &do_deltam_scaling_, &do_partials_, &do_planpar_, &do_enhanced_ps_, do_sources_up_.dataFirst(), do_sources_up_p_.dataFirst(), do_sources_dn_.dataFirst(), do_sources_dn_p_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), bb_input_.dataFirst(), &surfbb_, user_emissivity_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), mu1_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_up_.dataFirst(), wfine_up_.dataFirst(), hfine_up_.dataFirst(), xfine_dn_.dataFirst(), wfine_dn_.dataFirst(), hfine_dn_.dataFirst(), xfine_up_p_.dataFirst(), wfine_up_p_.dataFirst(), hfine_up_p_.dataFirst(), xfine_dn_p_.dataFirst(), wfine_dn_p_.dataFirst(), hfine_dn_p_.dataFirst(), intensity_dta_up_.dataFirst(), intensity_dts_.dataFirst(), intensity_dta_dn_.dataFirst(), cumsource_up_.dataFirst(), cumsource_dn_.dataFirst(), tcom1_.dataFirst(), lostrans_up_.dataFirst(), lostrans_up_p_.dataFirst());
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "Fo_Thermal_Rtcalcs_I:" << std::endl
      << "         maxgeoms: " << maxgeoms()  << std::endl
      << "        maxlayers: " << maxlayers()  << std::endl
      << "      maxpartials: " << maxpartials()  << std::endl
      << "          maxfine: " << maxfine()  << std::endl
      << "  max_user_levels: " << max_user_levels()  << std::endl
      << "      do_thermset: " << do_thermset()  << std::endl
      << "do_deltam_scaling: " << do_deltam_scaling()  << std::endl
      << "      do_partials: " << do_partials()  << std::endl
      << "       do_planpar: " << do_planpar()  << std::endl
      << "   do_enhanced_ps: " << do_enhanced_ps()  << std::endl
      << "    do_sources_dn: " << std::endl << do_sources_dn()  << std::endl
      << "  do_sources_dn_p: " << std::endl << do_sources_dn_p()  << std::endl
      << "           ngeoms: " << ngeoms()  << std::endl
      << "          nlayers: " << nlayers()  << std::endl
      << "        nfinedivs: " << std::endl << nfinedivs()  << std::endl
      << "    n_user_levels: " << n_user_levels()  << std::endl
      << "      user_levels: " << std::endl << user_levels()  << std::endl
      << "        npartials: " << npartials()  << std::endl
      << "      nfinedivs_p: " << std::endl << nfinedivs_p()  << std::endl
      << " partial_outindex: " << std::endl << partial_outindex()  << std::endl
      << "  partial_outflag: " << std::endl << partial_outflag()  << std::endl
      << " partial_layeridx: " << std::endl << partial_layeridx()  << std::endl
      << "         bb_input: " << std::endl << bb_input()  << std::endl
      << "       extinction: " << std::endl << extinction()  << std::endl
      << "          deltaus: " << std::endl << deltaus()  << std::endl
      << "            omega: " << std::endl << omega()  << std::endl
      << "         truncfac: " << std::endl << truncfac()  << std::endl
      << "              mu1: " << std::endl << mu1()  << std::endl
      << "       losw_paths: " << std::endl << losw_paths()  << std::endl
      << "       losp_paths: " << std::endl << losp_paths()  << std::endl
      << "            xfine: " << std::endl << xfine()  << std::endl
      << "            wfine: " << std::endl << wfine()  << std::endl
      << "            hfine: " << std::endl << hfine()  << std::endl
      << "          xfine_p: " << std::endl << xfine_p()  << std::endl
      << "          wfine_p: " << std::endl << wfine_p()  << std::endl
      << "          hfine_p: " << std::endl << hfine_p()  << std::endl
      << " intensity_dta_dn: " << std::endl << intensity_dta_dn()  << std::endl
      << "     cumsource_dn: " << std::endl << cumsource_dn()  << std::endl
      << "            tcom1: " << std::endl << tcom1()  << std::endl
      << "    do_sources_up: " << std::endl << do_sources_up()  << std::endl
      << "  do_sources_up_p: " << std::endl << do_sources_up_p()  << std::endl
      << "           surfbb: " << surfbb()  << std::endl
      << "  user_emissivity: " << std::endl << user_emissivity()  << std::endl
      << " intensity_dta_up: " << std::endl << intensity_dta_up()  << std::endl
      << "    intensity_dts: " << std::endl << intensity_dts()  << std::endl
      << "     cumsource_up: " << std::endl << cumsource_up()  << std::endl
      << "      lostrans_up: " << std::endl << lostrans_up()  << std::endl
      << "    lostrans_up_p: " << std::endl << lostrans_up_p()  << std::endl
      << "     do_upwelling: " << do_upwelling()  << std::endl
      << "     do_dnwelling: " << do_dnwelling()  << std::endl
      << "         xfine_up: " << std::endl << xfine_up()  << std::endl
      << "         wfine_up: " << std::endl << wfine_up()  << std::endl
      << "         hfine_up: " << std::endl << hfine_up()  << std::endl
      << "         xfine_dn: " << std::endl << xfine_dn()  << std::endl
      << "         wfine_dn: " << std::endl << wfine_dn()  << std::endl
      << "         hfine_dn: " << std::endl << hfine_dn()  << std::endl
      << "       xfine_up_p: " << std::endl << xfine_up_p()  << std::endl
      << "       wfine_up_p: " << std::endl << wfine_up_p()  << std::endl
      << "       hfine_up_p: " << std::endl << hfine_up_p()  << std::endl
      << "       xfine_dn_p: " << std::endl << xfine_dn_p()  << std::endl
      << "       wfine_dn_p: " << std::endl << wfine_dn_p()  << std::endl
      << "       hfine_dn_p: " << std::endl << hfine_dn_p()  << std::endl;

  }

private:
  int maxgeoms_;
  int maxlayers_;
  int maxpartials_;
  int maxfine_;
  int max_user_levels_;
  bool do_thermset_;
  bool do_deltam_scaling_;
  bool do_partials_;
  bool do_planpar_;
  bool do_enhanced_ps_;
  blitz::Array<bool, 2> do_sources_dn_;
  blitz::Array<bool, 2> do_sources_dn_p_;
  int ngeoms_;
  int nlayers_;
  blitz::Array<int, 2> nfinedivs_;
  int n_user_levels_;
  blitz::Array<int, 1> user_levels_;
  int npartials_;
  blitz::Array<int, 2> nfinedivs_p_;
  blitz::Array<int, 1> partial_outindex_;
  blitz::Array<bool, 1> partial_outflag_;
  blitz::Array<int, 1> partial_layeridx_;
  blitz::Array<double, 1> bb_input_;
  blitz::Array<double, 1> extinction_;
  blitz::Array<double, 1> deltaus_;
  blitz::Array<double, 1> omega_;
  blitz::Array<double, 1> truncfac_;
  blitz::Array<double, 1> mu1_;
  blitz::Array<double, 2> losw_paths_;
  blitz::Array<double, 2> losp_paths_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> hfine_;
  blitz::Array<double, 3> xfine_p_;
  blitz::Array<double, 3> wfine_p_;
  blitz::Array<double, 3> hfine_p_;
  blitz::Array<double, 2> intensity_dta_dn_;
  blitz::Array<double, 2> cumsource_dn_;
  blitz::Array<double, 2> tcom1_;
  blitz::Array<bool, 2> do_sources_up_;
  blitz::Array<bool, 2> do_sources_up_p_;
  double surfbb_;
  blitz::Array<double, 1> user_emissivity_;
  blitz::Array<double, 2> intensity_dta_up_;
  blitz::Array<double, 2> intensity_dts_;
  blitz::Array<double, 2> cumsource_up_;
  blitz::Array<double, 2> lostrans_up_;
  blitz::Array<double, 2> lostrans_up_p_;
  bool do_upwelling_;
  bool do_dnwelling_;
  blitz::Array<double, 3> xfine_up_;
  blitz::Array<double, 3> wfine_up_;
  blitz::Array<double, 3> hfine_up_;
  blitz::Array<double, 3> xfine_dn_;
  blitz::Array<double, 3> wfine_dn_;
  blitz::Array<double, 3> hfine_dn_;
  blitz::Array<double, 3> xfine_up_p_;
  blitz::Array<double, 3> wfine_up_p_;
  blitz::Array<double, 3> hfine_up_p_;
  blitz::Array<double, 3> xfine_dn_p_;
  blitz::Array<double, 3> wfine_dn_p_;
  blitz::Array<double, 3> hfine_dn_p_;

  // Serialization support
  Fo_Thermal_Rtcalcs_I() {}
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

//-----------------------------------------------------------------------
// Links to module: "fo_thermal_rtcalcs_ilps_m" in file: "FO_Thermal_RTCalcs_ILPS.f90"
//-----------------------------------------------------------------------

extern "C" {
  void fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_dn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* max_user_levels_in, const int* max_atmoswfs_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const bool* do_sources_dn_in, const bool* do_sources_dn_p_in, const bool* do_profilewfs_in, const bool* lvaryflags_in, const int* lvarynums_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* bb_input_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* l_extinction_in, const double* l_deltaus_in, const double* l_omega_in, const double* l_truncfac_in, const double* mu1_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_in, const double* wfine_in, const double* hfine_in, const double* xfine_p_in, const double* wfine_p_in, const double* hfine_p_in, const double* intensity_dta_dn_in, const double* lp_jacobians_dta_dn_in, const double* tcom1_in, const double* l_tcom1_in);
  void fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_up_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* max_user_levels_in, const int* max_atmoswfs_in, const int* max_surfacewfs_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const bool* do_sources_up_in, const bool* do_sources_up_p_in, const bool* do_profilewfs_in, const bool* do_surfacewfs_in, const bool* lvaryflags_in, const int* lvarynums_in, const int* n_surfacewfs_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* bb_input_in, const double* surfbb_in, const double* user_emissivity_in, const double* ls_user_emissivity_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* l_extinction_in, const double* l_deltaus_in, const double* l_omega_in, const double* l_truncfac_in, const double* mu1_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_in, const double* wfine_in, const double* hfine_in, const double* xfine_p_in, const double* wfine_p_in, const double* hfine_p_in, const double* intensity_dta_up_in, const double* intensity_dts_in, const double* lp_jacobians_dta_up_in, const double* lp_jacobians_dts_up_in, const double* ls_jacobians_dts_in, const double* tcom1_in, const double* l_tcom1_in, const double* lostrans_up_in, const double* lostrans_up_p_in, const double* l_lostrans_up_in, const double* l_lostrans_up_p_in);
  void fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_updn_wrap(const int* maxgeoms_in, const int* maxlayers_in, const int* maxpartials_in, const int* maxfine_in, const int* max_user_levels_in, const int* max_atmoswfs_in, const int* max_surfacewfs_in, const bool* do_upwelling_in, const bool* do_dnwelling_in, const bool* do_thermset_in, const bool* do_deltam_scaling_in, const bool* do_partials_in, const bool* do_planpar_in, const bool* do_enhanced_ps_in, const bool* do_sources_up_in, const bool* do_sources_up_p_in, const bool* do_sources_dn_in, const bool* do_sources_dn_p_in, const bool* do_profilewfs_in, const bool* do_surfacewfs_in, const bool* lvaryflags_in, const int* lvarynums_in, const int* n_surfacewfs_in, const int* ngeoms_in, const int* nlayers_in, const int* nfinedivs_in, const int* n_user_levels_in, const int* user_levels_in, const int* npartials_in, const int* nfinedivs_p_in, const int* partial_outindex_in, const bool* partial_outflag_in, const int* partial_layeridx_in, const double* bb_input_in, const double* surfbb_in, const double* user_emissivity_in, const double* ls_user_emissivity_in, const double* extinction_in, const double* deltaus_in, const double* omega_in, const double* truncfac_in, const double* l_extinction_in, const double* l_deltaus_in, const double* l_omega_in, const double* l_truncfac_in, const double* mu1_in, const double* losw_paths_in, const double* losp_paths_in, const double* xfine_up_in, const double* wfine_up_in, const double* hfine_up_in, const double* xfine_dn_in, const double* wfine_dn_in, const double* hfine_dn_in, const double* xfine_up_p_in, const double* wfine_up_p_in, const double* hfine_up_p_in, const double* xfine_dn_p_in, const double* wfine_dn_p_in, const double* hfine_dn_p_in, const double* intensity_dta_up_in, const double* intensity_dts_in, const double* intensity_dta_dn_in, const double* lp_jacobians_dta_up_in, const double* lp_jacobians_dts_up_in, const double* ls_jacobians_dts_in, const double* lp_jacobians_dta_dn_in, const double* tcom1_in, const double* l_tcom1_in, const double* lostrans_up_in, const double* lostrans_up_p_in, const double* l_lostrans_up_in, const double* l_lostrans_up_p_in);
}

class Fo_Thermal_Rtcalcs_Ilps : public Printable<Fo_Thermal_Rtcalcs_Ilps> {

public:
  Fo_Thermal_Rtcalcs_Ilps(const int& maxgeoms_in, const int& maxlayers_in, const int& maxpartials_in, const int& maxfine_in, const int& max_user_levels_in, const int& max_atmoswfs_in, const int& ngeoms_in, const int& nlayers_in, const int& n_user_levels_in, const int& npartials_in, const int& max_surfacewfs_in, const int& n_surfacewfs_in) : maxgeoms_(maxgeoms_in), maxlayers_(maxlayers_in), maxpartials_(maxpartials_in), maxfine_(maxfine_in), max_user_levels_(max_user_levels_in), max_atmoswfs_(max_atmoswfs_in), ngeoms_(ngeoms_in), nlayers_(nlayers_in), n_user_levels_(n_user_levels_in), npartials_(npartials_in), max_surfacewfs_(max_surfacewfs_in), n_surfacewfs_(n_surfacewfs_in) 
  { 
    do_thermset_ = false;
    do_deltam_scaling_ = false;
    do_partials_ = false;
    do_planpar_ = false;
    do_enhanced_ps_ = false;
    do_sources_dn_.reference( blitz::Array<bool, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_dn_ = false;
    do_sources_dn_p_.reference( blitz::Array<bool, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_dn_p_ = false;
    do_profilewfs_ = false;
    lvaryflags_.reference( blitz::Array<bool, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    lvaryflags_ = false;
    lvarynums_.reference( blitz::Array<int, 1>(maxlayers_, blitz::ColumnMajorArray<1>()) );
    lvarynums_ = 0;
    nfinedivs_.reference( blitz::Array<int, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_ = 0;
    user_levels_.reference( blitz::Array<int, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    user_levels_ = 0;
    nfinedivs_p_.reference( blitz::Array<int, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    nfinedivs_p_ = 0;
    partial_outindex_.reference( blitz::Array<int, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    partial_outindex_ = 0;
    partial_outflag_.reference( blitz::Array<bool, 1>(max_user_levels_, blitz::ColumnMajorArray<1>()) );
    partial_outflag_ = false;
    partial_layeridx_.reference( blitz::Array<int, 1>(maxpartials_, blitz::ColumnMajorArray<1>()) );
    partial_layeridx_ = 0;
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
    losw_paths_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losw_paths_ = 0;
    losp_paths_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    losp_paths_ = 0;
    xfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_ = 0;
    wfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_ = 0;
    hfine_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_ = 0;
    xfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_p_ = 0;
    wfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_p_ = 0;
    hfine_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_p_ = 0;
    intensity_dta_dn_.reference( blitz::Array<double, 2>(max_user_levels_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    intensity_dta_dn_ = 0;
    lp_jacobians_dta_dn_.reference( blitz::Array<double, 4>(max_user_levels_, maxgeoms_, maxlayers_, max_atmoswfs_, blitz::ColumnMajorArray<4>()) );
    lp_jacobians_dta_dn_ = 0;
    tcom1_.reference( blitz::Array<double, 2>(maxlayers_, 2, blitz::ColumnMajorArray<2>()) );
    tcom1_ = 0;
    l_tcom1_.reference( blitz::Array<double, 3>(maxlayers_, 2, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    l_tcom1_ = 0;
    do_sources_up_.reference( blitz::Array<bool, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_up_ = false;
    do_sources_up_p_.reference( blitz::Array<bool, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    do_sources_up_p_ = false;
    do_surfacewfs_ = false;
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
    lostrans_up_.reference( blitz::Array<double, 2>(maxlayers_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    lostrans_up_ = 0;
    lostrans_up_p_.reference( blitz::Array<double, 2>(maxpartials_, maxgeoms_, blitz::ColumnMajorArray<2>()) );
    lostrans_up_p_ = 0;
    l_lostrans_up_.reference( blitz::Array<double, 3>(maxlayers_, maxgeoms_, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    l_lostrans_up_ = 0;
    l_lostrans_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxgeoms_, max_atmoswfs_, blitz::ColumnMajorArray<3>()) );
    l_lostrans_up_p_ = 0;
    do_upwelling_ = false;
    do_dnwelling_ = false;
    xfine_up_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_up_ = 0;
    wfine_up_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_up_ = 0;
    hfine_up_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_up_ = 0;
    xfine_dn_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_dn_ = 0;
    wfine_dn_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_dn_ = 0;
    hfine_dn_.reference( blitz::Array<double, 3>(maxlayers_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_dn_ = 0;
    xfine_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_up_p_ = 0;
    wfine_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_up_p_ = 0;
    hfine_up_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_up_p_ = 0;
    xfine_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    xfine_dn_p_ = 0;
    wfine_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    wfine_dn_p_ = 0;
    hfine_dn_p_.reference( blitz::Array<double, 3>(maxpartials_, maxfine_, maxgeoms_, blitz::ColumnMajorArray<3>()) );
    hfine_dn_p_ = 0;
    // Initialize type pointers
    
  }

  virtual ~Fo_Thermal_Rtcalcs_Ilps() = default;

  const int& maxgeoms() const {
    return maxgeoms_;
  }

  

  const int& maxlayers() const {
    return maxlayers_;
  }

  

  const int& maxpartials() const {
    return maxpartials_;
  }

  

  const int& maxfine() const {
    return maxfine_;
  }

  

  const int& max_user_levels() const {
    return max_user_levels_;
  }

  

  const int& max_atmoswfs() const {
    return max_atmoswfs_;
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

  

  const bool& do_partials() const {
    return do_partials_;
  }

  void do_partials(const bool& do_partials_in) {
    do_partials_ = do_partials_in;
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

  

  const blitz::Array<bool, 2>& do_sources_dn() const {
    return do_sources_dn_;
  }

  void do_sources_dn(const blitz::Array<bool, 2>& do_sources_dn_in) {
    do_sources_dn_ = do_sources_dn_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_dn_p() const {
    return do_sources_dn_p_;
  }

  void do_sources_dn_p(const blitz::Array<bool, 2>& do_sources_dn_p_in) {
    do_sources_dn_p_ = do_sources_dn_p_in;
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

  

  const int& npartials() const {
    return npartials_;
  }

  

  const blitz::Array<int, 2>& nfinedivs_p() const {
    return nfinedivs_p_;
  }

  void nfinedivs_p(const blitz::Array<int, 2>& nfinedivs_p_in) {
    nfinedivs_p_ = nfinedivs_p_in;
  }

  

  const blitz::Array<int, 1>& partial_outindex() const {
    return partial_outindex_;
  }

  void partial_outindex(const blitz::Array<int, 1>& partial_outindex_in) {
    partial_outindex_ = partial_outindex_in;
  }

  

  const blitz::Array<bool, 1>& partial_outflag() const {
    return partial_outflag_;
  }

  void partial_outflag(const blitz::Array<bool, 1>& partial_outflag_in) {
    partial_outflag_ = partial_outflag_in;
  }

  

  const blitz::Array<int, 1>& partial_layeridx() const {
    return partial_layeridx_;
  }

  void partial_layeridx(const blitz::Array<int, 1>& partial_layeridx_in) {
    partial_layeridx_ = partial_layeridx_in;
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

  

  const blitz::Array<double, 2>& losw_paths() const {
    return losw_paths_;
  }

  void losw_paths(const blitz::Array<double, 2>& losw_paths_in) {
    losw_paths_ = losw_paths_in;
  }

  

  const blitz::Array<double, 2>& losp_paths() const {
    return losp_paths_;
  }

  void losp_paths(const blitz::Array<double, 2>& losp_paths_in) {
    losp_paths_ = losp_paths_in;
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

  

  const blitz::Array<double, 3>& hfine() const {
    return hfine_;
  }

  void hfine(const blitz::Array<double, 3>& hfine_in) {
    hfine_ = hfine_in;
  }

  

  const blitz::Array<double, 3>& xfine_p() const {
    return xfine_p_;
  }

  void xfine_p(const blitz::Array<double, 3>& xfine_p_in) {
    xfine_p_ = xfine_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_p() const {
    return wfine_p_;
  }

  void wfine_p(const blitz::Array<double, 3>& wfine_p_in) {
    wfine_p_ = wfine_p_in;
  }

  

  const blitz::Array<double, 3>& hfine_p() const {
    return hfine_p_;
  }

  void hfine_p(const blitz::Array<double, 3>& hfine_p_in) {
    hfine_p_ = hfine_p_in;
  }

  

  const blitz::Array<double, 2>& intensity_dta_dn() const {
    return intensity_dta_dn_;
  }

  

  const blitz::Array<double, 4>& lp_jacobians_dta_dn() const {
    return lp_jacobians_dta_dn_;
  }

  

  const blitz::Array<double, 2>& tcom1() const {
    return tcom1_;
  }

  void tcom1(const blitz::Array<double, 2>& tcom1_in) {
    tcom1_ = tcom1_in;
  }

  

  const blitz::Array<double, 3>& l_tcom1() const {
    return l_tcom1_;
  }

  void l_tcom1(const blitz::Array<double, 3>& l_tcom1_in) {
    l_tcom1_ = l_tcom1_in;
  }

  

  const int& max_surfacewfs() const {
    return max_surfacewfs_;
  }

  

  const blitz::Array<bool, 2>& do_sources_up() const {
    return do_sources_up_;
  }

  void do_sources_up(const blitz::Array<bool, 2>& do_sources_up_in) {
    do_sources_up_ = do_sources_up_in;
  }

  

  const blitz::Array<bool, 2>& do_sources_up_p() const {
    return do_sources_up_p_;
  }

  void do_sources_up_p(const blitz::Array<bool, 2>& do_sources_up_p_in) {
    do_sources_up_p_ = do_sources_up_p_in;
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

  

  const blitz::Array<double, 2>& lostrans_up() const {
    return lostrans_up_;
  }

  

  const blitz::Array<double, 2>& lostrans_up_p() const {
    return lostrans_up_p_;
  }

  

  const blitz::Array<double, 3>& l_lostrans_up() const {
    return l_lostrans_up_;
  }

  

  const blitz::Array<double, 3>& l_lostrans_up_p() const {
    return l_lostrans_up_p_;
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

  

  const blitz::Array<double, 3>& xfine_up() const {
    return xfine_up_;
  }

  void xfine_up(const blitz::Array<double, 3>& xfine_up_in) {
    xfine_up_ = xfine_up_in;
  }

  

  const blitz::Array<double, 3>& wfine_up() const {
    return wfine_up_;
  }

  void wfine_up(const blitz::Array<double, 3>& wfine_up_in) {
    wfine_up_ = wfine_up_in;
  }

  

  const blitz::Array<double, 3>& hfine_up() const {
    return hfine_up_;
  }

  void hfine_up(const blitz::Array<double, 3>& hfine_up_in) {
    hfine_up_ = hfine_up_in;
  }

  

  const blitz::Array<double, 3>& xfine_dn() const {
    return xfine_dn_;
  }

  void xfine_dn(const blitz::Array<double, 3>& xfine_dn_in) {
    xfine_dn_ = xfine_dn_in;
  }

  

  const blitz::Array<double, 3>& wfine_dn() const {
    return wfine_dn_;
  }

  void wfine_dn(const blitz::Array<double, 3>& wfine_dn_in) {
    wfine_dn_ = wfine_dn_in;
  }

  

  const blitz::Array<double, 3>& hfine_dn() const {
    return hfine_dn_;
  }

  void hfine_dn(const blitz::Array<double, 3>& hfine_dn_in) {
    hfine_dn_ = hfine_dn_in;
  }

  

  const blitz::Array<double, 3>& xfine_up_p() const {
    return xfine_up_p_;
  }

  void xfine_up_p(const blitz::Array<double, 3>& xfine_up_p_in) {
    xfine_up_p_ = xfine_up_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_up_p() const {
    return wfine_up_p_;
  }

  void wfine_up_p(const blitz::Array<double, 3>& wfine_up_p_in) {
    wfine_up_p_ = wfine_up_p_in;
  }

  

  const blitz::Array<double, 3>& hfine_up_p() const {
    return hfine_up_p_;
  }

  void hfine_up_p(const blitz::Array<double, 3>& hfine_up_p_in) {
    hfine_up_p_ = hfine_up_p_in;
  }

  

  const blitz::Array<double, 3>& xfine_dn_p() const {
    return xfine_dn_p_;
  }

  void xfine_dn_p(const blitz::Array<double, 3>& xfine_dn_p_in) {
    xfine_dn_p_ = xfine_dn_p_in;
  }

  

  const blitz::Array<double, 3>& wfine_dn_p() const {
    return wfine_dn_p_;
  }

  void wfine_dn_p(const blitz::Array<double, 3>& wfine_dn_p_in) {
    wfine_dn_p_ = wfine_dn_p_in;
  }

  

  const blitz::Array<double, 3>& hfine_dn_p() const {
    return hfine_dn_p_;
  }

  void hfine_dn_p(const blitz::Array<double, 3>& hfine_dn_p_in) {
    hfine_dn_p_ = hfine_dn_p_in;
  }

  

  
  void dte_integral_ilps_dn() {
    
    
    fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_dn_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &max_user_levels_, &max_atmoswfs_, &do_thermset_, &do_deltam_scaling_, &do_partials_, &do_planpar_, &do_enhanced_ps_, do_sources_dn_.dataFirst(), do_sources_dn_p_.dataFirst(), &do_profilewfs_, lvaryflags_.dataFirst(), lvarynums_.dataFirst(), &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), bb_input_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), l_extinction_.dataFirst(), l_deltaus_.dataFirst(), l_omega_.dataFirst(), l_truncfac_.dataFirst(), mu1_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), hfine_.dataFirst(), xfine_p_.dataFirst(), wfine_p_.dataFirst(), hfine_p_.dataFirst(), intensity_dta_dn_.dataFirst(), lp_jacobians_dta_dn_.dataFirst(), tcom1_.dataFirst(), l_tcom1_.dataFirst());
    
  }
void dte_integral_ilps_up() {
    
    
    fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_up_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &max_user_levels_, &max_atmoswfs_, &max_surfacewfs_, &do_thermset_, &do_deltam_scaling_, &do_partials_, &do_planpar_, &do_enhanced_ps_, do_sources_up_.dataFirst(), do_sources_up_p_.dataFirst(), &do_profilewfs_, &do_surfacewfs_, lvaryflags_.dataFirst(), lvarynums_.dataFirst(), &n_surfacewfs_, &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), bb_input_.dataFirst(), &surfbb_, user_emissivity_.dataFirst(), ls_user_emissivity_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), l_extinction_.dataFirst(), l_deltaus_.dataFirst(), l_omega_.dataFirst(), l_truncfac_.dataFirst(), mu1_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_.dataFirst(), wfine_.dataFirst(), hfine_.dataFirst(), xfine_p_.dataFirst(), wfine_p_.dataFirst(), hfine_p_.dataFirst(), intensity_dta_up_.dataFirst(), intensity_dts_.dataFirst(), lp_jacobians_dta_up_.dataFirst(), lp_jacobians_dts_up_.dataFirst(), ls_jacobians_dts_.dataFirst(), tcom1_.dataFirst(), l_tcom1_.dataFirst(), lostrans_up_.dataFirst(), lostrans_up_p_.dataFirst(), l_lostrans_up_.dataFirst(), l_lostrans_up_p_.dataFirst());
    
  }
void dte_integral_ilps_updn() {
    
    
    fo_thermal_rtcalcs_ilps_m_dte_integral_ilps_updn_wrap(&maxgeoms_, &maxlayers_, &maxpartials_, &maxfine_, &max_user_levels_, &max_atmoswfs_, &max_surfacewfs_, &do_upwelling_, &do_dnwelling_, &do_thermset_, &do_deltam_scaling_, &do_partials_, &do_planpar_, &do_enhanced_ps_, do_sources_up_.dataFirst(), do_sources_up_p_.dataFirst(), do_sources_dn_.dataFirst(), do_sources_dn_p_.dataFirst(), &do_profilewfs_, &do_surfacewfs_, lvaryflags_.dataFirst(), lvarynums_.dataFirst(), &n_surfacewfs_, &ngeoms_, &nlayers_, nfinedivs_.dataFirst(), &n_user_levels_, user_levels_.dataFirst(), &npartials_, nfinedivs_p_.dataFirst(), partial_outindex_.dataFirst(), partial_outflag_.dataFirst(), partial_layeridx_.dataFirst(), bb_input_.dataFirst(), &surfbb_, user_emissivity_.dataFirst(), ls_user_emissivity_.dataFirst(), extinction_.dataFirst(), deltaus_.dataFirst(), omega_.dataFirst(), truncfac_.dataFirst(), l_extinction_.dataFirst(), l_deltaus_.dataFirst(), l_omega_.dataFirst(), l_truncfac_.dataFirst(), mu1_.dataFirst(), losw_paths_.dataFirst(), losp_paths_.dataFirst(), xfine_up_.dataFirst(), wfine_up_.dataFirst(), hfine_up_.dataFirst(), xfine_dn_.dataFirst(), wfine_dn_.dataFirst(), hfine_dn_.dataFirst(), xfine_up_p_.dataFirst(), wfine_up_p_.dataFirst(), hfine_up_p_.dataFirst(), xfine_dn_p_.dataFirst(), wfine_dn_p_.dataFirst(), hfine_dn_p_.dataFirst(), intensity_dta_up_.dataFirst(), intensity_dts_.dataFirst(), intensity_dta_dn_.dataFirst(), lp_jacobians_dta_up_.dataFirst(), lp_jacobians_dts_up_.dataFirst(), ls_jacobians_dts_.dataFirst(), lp_jacobians_dta_dn_.dataFirst(), tcom1_.dataFirst(), l_tcom1_.dataFirst(), lostrans_up_.dataFirst(), lostrans_up_p_.dataFirst(), l_lostrans_up_.dataFirst(), l_lostrans_up_p_.dataFirst());
    
  }

  virtual void print(std::ostream &output_stream) const {
    output_stream << "Fo_Thermal_Rtcalcs_Ilps:" << std::endl
      << "           maxgeoms: " << maxgeoms()  << std::endl
      << "          maxlayers: " << maxlayers()  << std::endl
      << "        maxpartials: " << maxpartials()  << std::endl
      << "            maxfine: " << maxfine()  << std::endl
      << "    max_user_levels: " << max_user_levels()  << std::endl
      << "       max_atmoswfs: " << max_atmoswfs()  << std::endl
      << "        do_thermset: " << do_thermset()  << std::endl
      << "  do_deltam_scaling: " << do_deltam_scaling()  << std::endl
      << "        do_partials: " << do_partials()  << std::endl
      << "         do_planpar: " << do_planpar()  << std::endl
      << "     do_enhanced_ps: " << do_enhanced_ps()  << std::endl
      << "      do_sources_dn: " << std::endl << do_sources_dn()  << std::endl
      << "    do_sources_dn_p: " << std::endl << do_sources_dn_p()  << std::endl
      << "      do_profilewfs: " << do_profilewfs()  << std::endl
      << "         lvaryflags: " << std::endl << lvaryflags()  << std::endl
      << "          lvarynums: " << std::endl << lvarynums()  << std::endl
      << "             ngeoms: " << ngeoms()  << std::endl
      << "            nlayers: " << nlayers()  << std::endl
      << "          nfinedivs: " << std::endl << nfinedivs()  << std::endl
      << "      n_user_levels: " << n_user_levels()  << std::endl
      << "        user_levels: " << std::endl << user_levels()  << std::endl
      << "          npartials: " << npartials()  << std::endl
      << "        nfinedivs_p: " << std::endl << nfinedivs_p()  << std::endl
      << "   partial_outindex: " << std::endl << partial_outindex()  << std::endl
      << "    partial_outflag: " << std::endl << partial_outflag()  << std::endl
      << "   partial_layeridx: " << std::endl << partial_layeridx()  << std::endl
      << "           bb_input: " << std::endl << bb_input()  << std::endl
      << "         extinction: " << std::endl << extinction()  << std::endl
      << "            deltaus: " << std::endl << deltaus()  << std::endl
      << "              omega: " << std::endl << omega()  << std::endl
      << "           truncfac: " << std::endl << truncfac()  << std::endl
      << "       l_extinction: " << std::endl << l_extinction()  << std::endl
      << "          l_deltaus: " << std::endl << l_deltaus()  << std::endl
      << "            l_omega: " << std::endl << l_omega()  << std::endl
      << "         l_truncfac: " << std::endl << l_truncfac()  << std::endl
      << "                mu1: " << std::endl << mu1()  << std::endl
      << "         losw_paths: " << std::endl << losw_paths()  << std::endl
      << "         losp_paths: " << std::endl << losp_paths()  << std::endl
      << "              xfine: " << std::endl << xfine()  << std::endl
      << "              wfine: " << std::endl << wfine()  << std::endl
      << "              hfine: " << std::endl << hfine()  << std::endl
      << "            xfine_p: " << std::endl << xfine_p()  << std::endl
      << "            wfine_p: " << std::endl << wfine_p()  << std::endl
      << "            hfine_p: " << std::endl << hfine_p()  << std::endl
      << "   intensity_dta_dn: " << std::endl << intensity_dta_dn()  << std::endl
      << "lp_jacobians_dta_dn: " << std::endl << lp_jacobians_dta_dn()  << std::endl
      << "              tcom1: " << std::endl << tcom1()  << std::endl
      << "            l_tcom1: " << std::endl << l_tcom1()  << std::endl
      << "     max_surfacewfs: " << max_surfacewfs()  << std::endl
      << "      do_sources_up: " << std::endl << do_sources_up()  << std::endl
      << "    do_sources_up_p: " << std::endl << do_sources_up_p()  << std::endl
      << "      do_surfacewfs: " << do_surfacewfs()  << std::endl
      << "       n_surfacewfs: " << n_surfacewfs()  << std::endl
      << "             surfbb: " << surfbb()  << std::endl
      << "    user_emissivity: " << std::endl << user_emissivity()  << std::endl
      << " ls_user_emissivity: " << std::endl << ls_user_emissivity()  << std::endl
      << "   intensity_dta_up: " << std::endl << intensity_dta_up()  << std::endl
      << "      intensity_dts: " << std::endl << intensity_dts()  << std::endl
      << "lp_jacobians_dta_up: " << std::endl << lp_jacobians_dta_up()  << std::endl
      << "lp_jacobians_dts_up: " << std::endl << lp_jacobians_dts_up()  << std::endl
      << "   ls_jacobians_dts: " << std::endl << ls_jacobians_dts()  << std::endl
      << "        lostrans_up: " << std::endl << lostrans_up()  << std::endl
      << "      lostrans_up_p: " << std::endl << lostrans_up_p()  << std::endl
      << "      l_lostrans_up: " << std::endl << l_lostrans_up()  << std::endl
      << "    l_lostrans_up_p: " << std::endl << l_lostrans_up_p()  << std::endl
      << "       do_upwelling: " << do_upwelling()  << std::endl
      << "       do_dnwelling: " << do_dnwelling()  << std::endl
      << "           xfine_up: " << std::endl << xfine_up()  << std::endl
      << "           wfine_up: " << std::endl << wfine_up()  << std::endl
      << "           hfine_up: " << std::endl << hfine_up()  << std::endl
      << "           xfine_dn: " << std::endl << xfine_dn()  << std::endl
      << "           wfine_dn: " << std::endl << wfine_dn()  << std::endl
      << "           hfine_dn: " << std::endl << hfine_dn()  << std::endl
      << "         xfine_up_p: " << std::endl << xfine_up_p()  << std::endl
      << "         wfine_up_p: " << std::endl << wfine_up_p()  << std::endl
      << "         hfine_up_p: " << std::endl << hfine_up_p()  << std::endl
      << "         xfine_dn_p: " << std::endl << xfine_dn_p()  << std::endl
      << "         wfine_dn_p: " << std::endl << wfine_dn_p()  << std::endl
      << "         hfine_dn_p: " << std::endl << hfine_dn_p()  << std::endl;

  }

private:
  int maxgeoms_;
  int maxlayers_;
  int maxpartials_;
  int maxfine_;
  int max_user_levels_;
  int max_atmoswfs_;
  bool do_thermset_;
  bool do_deltam_scaling_;
  bool do_partials_;
  bool do_planpar_;
  bool do_enhanced_ps_;
  blitz::Array<bool, 2> do_sources_dn_;
  blitz::Array<bool, 2> do_sources_dn_p_;
  bool do_profilewfs_;
  blitz::Array<bool, 1> lvaryflags_;
  blitz::Array<int, 1> lvarynums_;
  int ngeoms_;
  int nlayers_;
  blitz::Array<int, 2> nfinedivs_;
  int n_user_levels_;
  blitz::Array<int, 1> user_levels_;
  int npartials_;
  blitz::Array<int, 2> nfinedivs_p_;
  blitz::Array<int, 1> partial_outindex_;
  blitz::Array<bool, 1> partial_outflag_;
  blitz::Array<int, 1> partial_layeridx_;
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
  blitz::Array<double, 2> losw_paths_;
  blitz::Array<double, 2> losp_paths_;
  blitz::Array<double, 3> xfine_;
  blitz::Array<double, 3> wfine_;
  blitz::Array<double, 3> hfine_;
  blitz::Array<double, 3> xfine_p_;
  blitz::Array<double, 3> wfine_p_;
  blitz::Array<double, 3> hfine_p_;
  blitz::Array<double, 2> intensity_dta_dn_;
  blitz::Array<double, 4> lp_jacobians_dta_dn_;
  blitz::Array<double, 2> tcom1_;
  blitz::Array<double, 3> l_tcom1_;
  int max_surfacewfs_;
  blitz::Array<bool, 2> do_sources_up_;
  blitz::Array<bool, 2> do_sources_up_p_;
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
  blitz::Array<double, 2> lostrans_up_;
  blitz::Array<double, 2> lostrans_up_p_;
  blitz::Array<double, 3> l_lostrans_up_;
  blitz::Array<double, 3> l_lostrans_up_p_;
  bool do_upwelling_;
  bool do_dnwelling_;
  blitz::Array<double, 3> xfine_up_;
  blitz::Array<double, 3> wfine_up_;
  blitz::Array<double, 3> hfine_up_;
  blitz::Array<double, 3> xfine_dn_;
  blitz::Array<double, 3> wfine_dn_;
  blitz::Array<double, 3> hfine_dn_;
  blitz::Array<double, 3> xfine_up_p_;
  blitz::Array<double, 3> wfine_up_p_;
  blitz::Array<double, 3> hfine_up_p_;
  blitz::Array<double, 3> xfine_dn_p_;
  blitz::Array<double, 3> wfine_dn_p_;
  blitz::Array<double, 3> hfine_dn_p_;

  // Serialization support
  Fo_Thermal_Rtcalcs_Ilps() {}
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);
  template<class Archive> void save(Archive & ar, const unsigned int version) const;
  template<class Archive> void load(Archive & ar, const unsigned int version);
};

}

// Serialization support
FP_EXPORT_KEY(Fo_Dtwpgeometry_Master)
FP_EXPORT_KEY(Fo_Sswpgeometry_Master)
FP_EXPORT_KEY(Fo_Scalarss_Rtcalcs_I)
FP_EXPORT_KEY(Fo_Scalarss_Rtcalcs_Ilps)
FP_EXPORT_KEY(Fo_Scalarss_Spherfuncs)
FP_EXPORT_KEY(Fo_Thermal_Rtcalcs_I)
FP_EXPORT_KEY(Fo_Thermal_Rtcalcs_Ilps)

#endif
