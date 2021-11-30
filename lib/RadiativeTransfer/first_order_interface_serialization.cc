#include "first_order_interface.h"
#include "fp_serialize_support.h"
#include "linear_algebra.h"

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION

////////
// Serialization functions of Fo_Dtwpgeometry_Master
template<class Archive>
void Fo_Dtwpgeometry_Master::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Fo_Dtwpgeometry_Master);
    ar 
       & FP_NVP_(maxgeoms)
       & FP_NVP_(maxlayers)
       & FP_NVP_(maxpartials)
       & FP_NVP_(maxfine)
       & FP_NVP_(dtr)
       & FP_NVP_(eradius)
       & FP_NVP_(do_upwelling)
       & FP_NVP_(do_planpar)
       & FP_NVP_(do_enhanced_ps)
       & FP_NVP_(do_partials)
       & FP_NVP_(ngeoms)
       & FP_NVP_(nlayers)
       & FP_NVP_(npartials)
       & FP_NVP_(nfine)
       & FP_NVP_(partial_layeridx)
       & FP_NVP_(heights)
       & FP_NVP_(alpha_boa)
       & FP_NVP_(partial_heights)
       & FP_NVP_(mu1)
       & FP_NVP_(radii)
       & FP_NVP_(losw_paths)
       & FP_NVP_(alpha)
       & FP_NVP_(sina)
       & FP_NVP_(cosa)
       & FP_NVP_(radii_p)
       & FP_NVP_(losp_paths)
       & FP_NVP_(alpha_p)
       & FP_NVP_(sina_p)
       & FP_NVP_(cosa_p)
       & FP_NVP_(nfinedivs)
       & FP_NVP_(xfine)
       & FP_NVP_(wfine)
       & FP_NVP_(radiifine)
       & FP_NVP_(alphafine)
       & FP_NVP_(sinfine)
       & FP_NVP_(cosfine)
       & FP_NVP_(nfinedivs_p)
       & FP_NVP_(xfine_p)
       & FP_NVP_(wfine_p)
       & FP_NVP_(radiifine_p)
       & FP_NVP_(alphafine_p)
       & FP_NVP_(sinfine_p)
       & FP_NVP_(cosfine_p)
       & FP_NVP_(fail)
       & FP_NVP_(message)
       & FP_NVP_(trace);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Fo_Dtwpgeometry_Master::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void Fo_Dtwpgeometry_Master::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  losw_paths_.reference(to_fortran(losw_paths_));
  alpha_.reference(to_fortran(alpha_));
  sina_.reference(to_fortran(sina_));
  cosa_.reference(to_fortran(cosa_));
  losp_paths_.reference(to_fortran(losp_paths_));
  alpha_p_.reference(to_fortran(alpha_p_));
  sina_p_.reference(to_fortran(sina_p_));
  cosa_p_.reference(to_fortran(cosa_p_));
  nfinedivs_.reference(to_fortran(nfinedivs_));
  xfine_.reference(to_fortran(xfine_));
  wfine_.reference(to_fortran(wfine_));
  radiifine_.reference(to_fortran(radiifine_));
  alphafine_.reference(to_fortran(alphafine_));
  sinfine_.reference(to_fortran(sinfine_));
  cosfine_.reference(to_fortran(cosfine_));
  nfinedivs_p_.reference(to_fortran(nfinedivs_p_));
  xfine_p_.reference(to_fortran(xfine_p_));
  wfine_p_.reference(to_fortran(wfine_p_));
  radiifine_p_.reference(to_fortran(radiifine_p_));
  alphafine_p_.reference(to_fortran(alphafine_p_));
  sinfine_p_.reference(to_fortran(sinfine_p_));
  cosfine_p_.reference(to_fortran(cosfine_p_));
}

////////
// Serialization functions of Fo_Sswpgeometry_Master
template<class Archive>
void Fo_Sswpgeometry_Master::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Fo_Sswpgeometry_Master);
    ar 
       & FP_NVP_(maxgeoms)
       & FP_NVP_(maxszas)
       & FP_NVP_(maxvzas)
       & FP_NVP_(maxazms)
       & FP_NVP_(maxlayers)
       & FP_NVP_(maxpartials)
       & FP_NVP_(maxfine)
       & FP_NVP_(do_obsgeom)
       & FP_NVP_(do_doublet)
       & FP_NVP_(do_chapman)
       & FP_NVP_(do_planpar)
       & FP_NVP_(do_enhanced_ps)
       & FP_NVP_(do_partials)
       & FP_NVP_(ngeoms)
       & FP_NVP_(nszas)
       & FP_NVP_(nvzas)
       & FP_NVP_(nazms)
       & FP_NVP_(nlayers)
       & FP_NVP_(nfine)
       & FP_NVP_(npartials)
       & FP_NVP_(partial_layeridx)
       & FP_NVP_(dtr)
       & FP_NVP_(pie)
       & FP_NVP_(vsign)
       & FP_NVP_(eradius)
       & FP_NVP_(nv_offset)
       & FP_NVP_(na_offset)
       & FP_NVP_(nd_offset)
       & FP_NVP_(heights)
       & FP_NVP_(partial_heights)
       & FP_NVP_(obsgeom_boa)
       & FP_NVP_(alpha_boa)
       & FP_NVP_(theta_boa)
       & FP_NVP_(phi_boa)
       & FP_NVP_(donadir)
       & FP_NVP_(raycon)
       & FP_NVP_(mu0)
       & FP_NVP_(mu1)
       & FP_NVP_(cosscat)
       & FP_NVP_(radii)
       & FP_NVP_(losw_paths)
       & FP_NVP_(alpha)
       & FP_NVP_(sina)
       & FP_NVP_(cosa)
       & FP_NVP_(sunpaths)
       & FP_NVP_(ntraverse)
       & FP_NVP_(chapfacs)
       & FP_NVP_(theta_all)
       & FP_NVP_(radii_p)
       & FP_NVP_(losp_paths)
       & FP_NVP_(alpha_p)
       & FP_NVP_(sina_p)
       & FP_NVP_(cosa_p)
       & FP_NVP_(sunpaths_p)
       & FP_NVP_(ntraverse_p)
       & FP_NVP_(chapfacs_p)
       & FP_NVP_(nfinedivs)
       & FP_NVP_(xfine)
       & FP_NVP_(wfine)
       & FP_NVP_(radiifine)
       & FP_NVP_(alphafine)
       & FP_NVP_(sinfine)
       & FP_NVP_(cosfine)
       & FP_NVP_(sunpathsfine)
       & FP_NVP_(ntraversefine)
       & FP_NVP_(nfinedivs_p)
       & FP_NVP_(xfine_p)
       & FP_NVP_(wfine_p)
       & FP_NVP_(radiifine_p)
       & FP_NVP_(alphafine_p)
       & FP_NVP_(sinfine_p)
       & FP_NVP_(cosfine_p)
       & FP_NVP_(sunpathsfine_p)
       & FP_NVP_(ntraversefine_p)
       & FP_NVP_(fail)
       & FP_NVP_(message)
       & FP_NVP_(trace);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Fo_Sswpgeometry_Master::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void Fo_Sswpgeometry_Master::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  na_offset_.reference(to_fortran(na_offset_));
  obsgeom_boa_.reference(to_fortran(obsgeom_boa_));
  losw_paths_.reference(to_fortran(losw_paths_));
  alpha_.reference(to_fortran(alpha_));
  sina_.reference(to_fortran(sina_));
  cosa_.reference(to_fortran(cosa_));
  sunpaths_.reference(to_fortran(sunpaths_));
  ntraverse_.reference(to_fortran(ntraverse_));
  chapfacs_.reference(to_fortran(chapfacs_));
  theta_all_.reference(to_fortran(theta_all_));
  losp_paths_.reference(to_fortran(losp_paths_));
  alpha_p_.reference(to_fortran(alpha_p_));
  sina_p_.reference(to_fortran(sina_p_));
  cosa_p_.reference(to_fortran(cosa_p_));
  sunpaths_p_.reference(to_fortran(sunpaths_p_));
  ntraverse_p_.reference(to_fortran(ntraverse_p_));
  chapfacs_p_.reference(to_fortran(chapfacs_p_));
  nfinedivs_.reference(to_fortran(nfinedivs_));
  xfine_.reference(to_fortran(xfine_));
  wfine_.reference(to_fortran(wfine_));
  radiifine_.reference(to_fortran(radiifine_));
  alphafine_.reference(to_fortran(alphafine_));
  sinfine_.reference(to_fortran(sinfine_));
  cosfine_.reference(to_fortran(cosfine_));
  sunpathsfine_.reference(to_fortran(sunpathsfine_));
  ntraversefine_.reference(to_fortran(ntraversefine_));
  nfinedivs_p_.reference(to_fortran(nfinedivs_p_));
  xfine_p_.reference(to_fortran(xfine_p_));
  wfine_p_.reference(to_fortran(wfine_p_));
  radiifine_p_.reference(to_fortran(radiifine_p_));
  alphafine_p_.reference(to_fortran(alphafine_p_));
  sinfine_p_.reference(to_fortran(sinfine_p_));
  cosfine_p_.reference(to_fortran(cosfine_p_));
  sunpathsfine_p_.reference(to_fortran(sunpathsfine_p_));
  ntraversefine_p_.reference(to_fortran(ntraversefine_p_));
}

////////
// Serialization functions of Fo_Scalarss_Rtcalcs_I
template<class Archive>
void Fo_Scalarss_Rtcalcs_I::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Fo_Scalarss_Rtcalcs_I);
    ar 
       & FP_NVP_(maxgeoms)
       & FP_NVP_(maxlayers)
       & FP_NVP_(maxpartials)
       & FP_NVP_(maxfine)
       & FP_NVP_(maxmoments_input)
       & FP_NVP_(max_user_levels)
       & FP_NVP_(do_deltam_scaling)
       & FP_NVP_(do_phasfunc)
       & FP_NVP_(do_partials)
       & FP_NVP_(do_planpar)
       & FP_NVP_(do_enhanced_ps)
       & FP_NVP_(flux)
       & FP_NVP_(do_sources_dn)
       & FP_NVP_(do_sources_dn_p)
       & FP_NVP_(ngeoms)
       & FP_NVP_(nlayers)
       & FP_NVP_(nfinedivs)
       & FP_NVP_(nmoments_input)
       & FP_NVP_(n_user_levels)
       & FP_NVP_(user_levels)
       & FP_NVP_(npartials)
       & FP_NVP_(nfinedivs_p)
       & FP_NVP_(partial_outindex)
       & FP_NVP_(partial_outflag)
       & FP_NVP_(partial_layeridx)
       & FP_NVP_(extinction)
       & FP_NVP_(deltaus)
       & FP_NVP_(omega)
       & FP_NVP_(truncfac)
       & FP_NVP_(phasmoms)
       & FP_NVP_(phasfunc_dn)
       & FP_NVP_(mu1)
       & FP_NVP_(legpoly_dn)
       & FP_NVP_(losw_paths)
       & FP_NVP_(losp_paths)
       & FP_NVP_(xfine)
       & FP_NVP_(wfine)
       & FP_NVP_(sunpaths)
       & FP_NVP_(ntraverse)
       & FP_NVP_(sunpathsfine)
       & FP_NVP_(ntraversefine)
       & FP_NVP_(xfine_p)
       & FP_NVP_(wfine_p)
       & FP_NVP_(sunpaths_p)
       & FP_NVP_(ntraverse_p)
       & FP_NVP_(sunpathsfine_p)
       & FP_NVP_(ntraversefine_p)
       & FP_NVP_(intensity_dn)
       & FP_NVP_(cumsource_dn)
       & FP_NVP_(lostrans_dn)
       & FP_NVP_(do_surface_leaving)
       & FP_NVP_(do_water_leaving)
       & FP_NVP_(do_sources_up)
       & FP_NVP_(do_sources_up_p)
       & FP_NVP_(phasfunc_up)
       & FP_NVP_(reflec)
       & FP_NVP_(slterm)
       & FP_NVP_(mu0)
       & FP_NVP_(legpoly_up)
       & FP_NVP_(intensity_up)
       & FP_NVP_(intensity_db)
       & FP_NVP_(cumsource_up)
       & FP_NVP_(cumtrans)
       & FP_NVP_(lostrans_up)
       & FP_NVP_(do_upwelling)
       & FP_NVP_(do_dnwelling)
       & FP_NVP_(xfine_up)
       & FP_NVP_(wfine_up)
       & FP_NVP_(sunpaths_up)
       & FP_NVP_(ntraverse_up)
       & FP_NVP_(sunpathsfine_up)
       & FP_NVP_(ntraversefine_up)
       & FP_NVP_(xfine_dn)
       & FP_NVP_(wfine_dn)
       & FP_NVP_(sunpaths_dn)
       & FP_NVP_(ntraverse_dn)
       & FP_NVP_(sunpathsfine_dn)
       & FP_NVP_(ntraversefine_dn)
       & FP_NVP_(xfine_up_p)
       & FP_NVP_(wfine_up_p)
       & FP_NVP_(sunpaths_up_p)
       & FP_NVP_(ntraverse_up_p)
       & FP_NVP_(sunpathsfine_up_p)
       & FP_NVP_(ntraversefine_up_p)
       & FP_NVP_(xfine_dn_p)
       & FP_NVP_(wfine_dn_p)
       & FP_NVP_(sunpaths_dn_p)
       & FP_NVP_(ntraverse_dn_p)
       & FP_NVP_(sunpathsfine_dn_p)
       & FP_NVP_(ntraversefine_dn_p);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Fo_Scalarss_Rtcalcs_I::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void Fo_Scalarss_Rtcalcs_I::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  do_sources_dn_.reference(to_fortran(do_sources_dn_));
  do_sources_dn_p_.reference(to_fortran(do_sources_dn_p_));
  nfinedivs_.reference(to_fortran(nfinedivs_));
  nfinedivs_p_.reference(to_fortran(nfinedivs_p_));
  phasmoms_.reference(to_fortran(phasmoms_));
  phasfunc_dn_.reference(to_fortran(phasfunc_dn_));
  legpoly_dn_.reference(to_fortran(legpoly_dn_));
  losw_paths_.reference(to_fortran(losw_paths_));
  losp_paths_.reference(to_fortran(losp_paths_));
  xfine_.reference(to_fortran(xfine_));
  wfine_.reference(to_fortran(wfine_));
  sunpaths_.reference(to_fortran(sunpaths_));
  ntraverse_.reference(to_fortran(ntraverse_));
  sunpathsfine_.reference(to_fortran(sunpathsfine_));
  ntraversefine_.reference(to_fortran(ntraversefine_));
  xfine_p_.reference(to_fortran(xfine_p_));
  wfine_p_.reference(to_fortran(wfine_p_));
  sunpaths_p_.reference(to_fortran(sunpaths_p_));
  ntraverse_p_.reference(to_fortran(ntraverse_p_));
  sunpathsfine_p_.reference(to_fortran(sunpathsfine_p_));
  ntraversefine_p_.reference(to_fortran(ntraversefine_p_));
  intensity_dn_.reference(to_fortran(intensity_dn_));
  cumsource_dn_.reference(to_fortran(cumsource_dn_));
  lostrans_dn_.reference(to_fortran(lostrans_dn_));
  do_sources_up_.reference(to_fortran(do_sources_up_));
  do_sources_up_p_.reference(to_fortran(do_sources_up_p_));
  phasfunc_up_.reference(to_fortran(phasfunc_up_));
  legpoly_up_.reference(to_fortran(legpoly_up_));
  intensity_up_.reference(to_fortran(intensity_up_));
  intensity_db_.reference(to_fortran(intensity_db_));
  cumsource_up_.reference(to_fortran(cumsource_up_));
  cumtrans_.reference(to_fortran(cumtrans_));
  lostrans_up_.reference(to_fortran(lostrans_up_));
  xfine_up_.reference(to_fortran(xfine_up_));
  wfine_up_.reference(to_fortran(wfine_up_));
  sunpaths_up_.reference(to_fortran(sunpaths_up_));
  ntraverse_up_.reference(to_fortran(ntraverse_up_));
  sunpathsfine_up_.reference(to_fortran(sunpathsfine_up_));
  ntraversefine_up_.reference(to_fortran(ntraversefine_up_));
  xfine_dn_.reference(to_fortran(xfine_dn_));
  wfine_dn_.reference(to_fortran(wfine_dn_));
  sunpaths_dn_.reference(to_fortran(sunpaths_dn_));
  ntraverse_dn_.reference(to_fortran(ntraverse_dn_));
  sunpathsfine_dn_.reference(to_fortran(sunpathsfine_dn_));
  ntraversefine_dn_.reference(to_fortran(ntraversefine_dn_));
  xfine_up_p_.reference(to_fortran(xfine_up_p_));
  wfine_up_p_.reference(to_fortran(wfine_up_p_));
  sunpaths_up_p_.reference(to_fortran(sunpaths_up_p_));
  ntraverse_up_p_.reference(to_fortran(ntraverse_up_p_));
  sunpathsfine_up_p_.reference(to_fortran(sunpathsfine_up_p_));
  ntraversefine_up_p_.reference(to_fortran(ntraversefine_up_p_));
  xfine_dn_p_.reference(to_fortran(xfine_dn_p_));
  wfine_dn_p_.reference(to_fortran(wfine_dn_p_));
  sunpaths_dn_p_.reference(to_fortran(sunpaths_dn_p_));
  ntraverse_dn_p_.reference(to_fortran(ntraverse_dn_p_));
  sunpathsfine_dn_p_.reference(to_fortran(sunpathsfine_dn_p_));
  ntraversefine_dn_p_.reference(to_fortran(ntraversefine_dn_p_));
}

////////
// Serialization functions of Fo_Scalarss_Rtcalcs_Ilps
template<class Archive>
void Fo_Scalarss_Rtcalcs_Ilps::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Fo_Scalarss_Rtcalcs_Ilps);
    ar 
       & FP_NVP_(maxgeoms)
       & FP_NVP_(maxlayers)
       & FP_NVP_(maxpartials)
       & FP_NVP_(maxfine)
       & FP_NVP_(maxmoments_input)
       & FP_NVP_(max_user_levels)
       & FP_NVP_(max_atmoswfs)
       & FP_NVP_(do_deltam_scaling)
       & FP_NVP_(do_phasfunc)
       & FP_NVP_(do_partials)
       & FP_NVP_(do_planpar)
       & FP_NVP_(do_enhanced_ps)
       & FP_NVP_(flux)
       & FP_NVP_(do_sources_dn)
       & FP_NVP_(do_sources_dn_p)
       & FP_NVP_(do_profilewfs)
       & FP_NVP_(lvaryflags)
       & FP_NVP_(lvarynums)
       & FP_NVP_(lvarymoms)
       & FP_NVP_(ngeoms)
       & FP_NVP_(nlayers)
       & FP_NVP_(nfinedivs)
       & FP_NVP_(nmoments_input)
       & FP_NVP_(n_user_levels)
       & FP_NVP_(user_levels)
       & FP_NVP_(npartials)
       & FP_NVP_(nfinedivs_p)
       & FP_NVP_(partial_outindex)
       & FP_NVP_(partial_outflag)
       & FP_NVP_(partial_layeridx)
       & FP_NVP_(extinction)
       & FP_NVP_(deltaus)
       & FP_NVP_(omega)
       & FP_NVP_(truncfac)
       & FP_NVP_(phasmoms)
       & FP_NVP_(phasfunc_dn)
       & FP_NVP_(l_extinction)
       & FP_NVP_(l_deltaus)
       & FP_NVP_(l_omega)
       & FP_NVP_(l_truncfac)
       & FP_NVP_(l_phasmoms)
       & FP_NVP_(l_phasfunc_dn)
       & FP_NVP_(mu1)
       & FP_NVP_(legpoly_dn)
       & FP_NVP_(losw_paths)
       & FP_NVP_(losp_paths)
       & FP_NVP_(xfine)
       & FP_NVP_(wfine)
       & FP_NVP_(sunpaths)
       & FP_NVP_(ntraverse)
       & FP_NVP_(sunpathsfine)
       & FP_NVP_(ntraversefine)
       & FP_NVP_(xfine_p)
       & FP_NVP_(wfine_p)
       & FP_NVP_(sunpaths_p)
       & FP_NVP_(ntraverse_p)
       & FP_NVP_(sunpathsfine_p)
       & FP_NVP_(ntraversefine_p)
       & FP_NVP_(intensity_dn)
       & FP_NVP_(lp_jacobians_dn)
       & FP_NVP_(lostrans_dn)
       & FP_NVP_(lp_lostrans_dn)
       & FP_NVP_(max_surfacewfs)
       & FP_NVP_(max_sleavewfs)
       & FP_NVP_(do_surface_leaving)
       & FP_NVP_(do_water_leaving)
       & FP_NVP_(do_sources_up)
       & FP_NVP_(do_sources_up_p)
       & FP_NVP_(do_surfacewfs)
       & FP_NVP_(do_sleavewfs)
       & FP_NVP_(n_reflecwfs)
       & FP_NVP_(n_sleavewfs)
       & FP_NVP_(n_surfacewfs)
       & FP_NVP_(phasfunc_up)
       & FP_NVP_(reflec)
       & FP_NVP_(slterm)
       & FP_NVP_(l_phasfunc_up)
       & FP_NVP_(ls_reflec)
       & FP_NVP_(lssl_slterm)
       & FP_NVP_(mu0)
       & FP_NVP_(legpoly_up)
       & FP_NVP_(intensity_up)
       & FP_NVP_(intensity_db)
       & FP_NVP_(lp_jacobians_up)
       & FP_NVP_(lp_jacobians_db)
       & FP_NVP_(ls_jacobians_db)
       & FP_NVP_(cumtrans)
       & FP_NVP_(lostrans_up)
       & FP_NVP_(lp_cumtrans)
       & FP_NVP_(lp_lostrans_up)
       & FP_NVP_(do_upwelling)
       & FP_NVP_(do_dnwelling)
       & FP_NVP_(xfine_up)
       & FP_NVP_(wfine_up)
       & FP_NVP_(sunpaths_up)
       & FP_NVP_(ntraverse_up)
       & FP_NVP_(sunpathsfine_up)
       & FP_NVP_(ntraversefine_up)
       & FP_NVP_(xfine_dn)
       & FP_NVP_(wfine_dn)
       & FP_NVP_(sunpaths_dn)
       & FP_NVP_(ntraverse_dn)
       & FP_NVP_(sunpathsfine_dn)
       & FP_NVP_(ntraversefine_dn)
       & FP_NVP_(xfine_up_p)
       & FP_NVP_(wfine_up_p)
       & FP_NVP_(sunpaths_up_p)
       & FP_NVP_(ntraverse_up_p)
       & FP_NVP_(sunpathsfine_up_p)
       & FP_NVP_(ntraversefine_up_p)
       & FP_NVP_(xfine_dn_p)
       & FP_NVP_(wfine_dn_p)
       & FP_NVP_(sunpaths_dn_p)
       & FP_NVP_(ntraverse_dn_p)
       & FP_NVP_(sunpathsfine_dn_p)
       & FP_NVP_(ntraversefine_dn_p);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Fo_Scalarss_Rtcalcs_Ilps::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void Fo_Scalarss_Rtcalcs_Ilps::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  do_sources_dn_.reference(to_fortran(do_sources_dn_));
  do_sources_dn_p_.reference(to_fortran(do_sources_dn_p_));
  lvarymoms_.reference(to_fortran(lvarymoms_));
  nfinedivs_.reference(to_fortran(nfinedivs_));
  nfinedivs_p_.reference(to_fortran(nfinedivs_p_));
  phasmoms_.reference(to_fortran(phasmoms_));
  phasfunc_dn_.reference(to_fortran(phasfunc_dn_));
  l_extinction_.reference(to_fortran(l_extinction_));
  l_deltaus_.reference(to_fortran(l_deltaus_));
  l_omega_.reference(to_fortran(l_omega_));
  l_truncfac_.reference(to_fortran(l_truncfac_));
  l_phasmoms_.reference(to_fortran(l_phasmoms_));
  l_phasfunc_dn_.reference(to_fortran(l_phasfunc_dn_));
  legpoly_dn_.reference(to_fortran(legpoly_dn_));
  losw_paths_.reference(to_fortran(losw_paths_));
  losp_paths_.reference(to_fortran(losp_paths_));
  xfine_.reference(to_fortran(xfine_));
  wfine_.reference(to_fortran(wfine_));
  sunpaths_.reference(to_fortran(sunpaths_));
  ntraverse_.reference(to_fortran(ntraverse_));
  sunpathsfine_.reference(to_fortran(sunpathsfine_));
  ntraversefine_.reference(to_fortran(ntraversefine_));
  xfine_p_.reference(to_fortran(xfine_p_));
  wfine_p_.reference(to_fortran(wfine_p_));
  sunpaths_p_.reference(to_fortran(sunpaths_p_));
  ntraverse_p_.reference(to_fortran(ntraverse_p_));
  sunpathsfine_p_.reference(to_fortran(sunpathsfine_p_));
  ntraversefine_p_.reference(to_fortran(ntraversefine_p_));
  intensity_dn_.reference(to_fortran(intensity_dn_));
  lp_jacobians_dn_.reference(to_fortran(lp_jacobians_dn_));
  lostrans_dn_.reference(to_fortran(lostrans_dn_));
  lp_lostrans_dn_.reference(to_fortran(lp_lostrans_dn_));
  do_sources_up_.reference(to_fortran(do_sources_up_));
  do_sources_up_p_.reference(to_fortran(do_sources_up_p_));
  phasfunc_up_.reference(to_fortran(phasfunc_up_));
  l_phasfunc_up_.reference(to_fortran(l_phasfunc_up_));
  ls_reflec_.reference(to_fortran(ls_reflec_));
  lssl_slterm_.reference(to_fortran(lssl_slterm_));
  legpoly_up_.reference(to_fortran(legpoly_up_));
  intensity_up_.reference(to_fortran(intensity_up_));
  intensity_db_.reference(to_fortran(intensity_db_));
  lp_jacobians_up_.reference(to_fortran(lp_jacobians_up_));
  lp_jacobians_db_.reference(to_fortran(lp_jacobians_db_));
  ls_jacobians_db_.reference(to_fortran(ls_jacobians_db_));
  cumtrans_.reference(to_fortran(cumtrans_));
  lostrans_up_.reference(to_fortran(lostrans_up_));
  lp_cumtrans_.reference(to_fortran(lp_cumtrans_));
  lp_lostrans_up_.reference(to_fortran(lp_lostrans_up_));
  xfine_up_.reference(to_fortran(xfine_up_));
  wfine_up_.reference(to_fortran(wfine_up_));
  sunpaths_up_.reference(to_fortran(sunpaths_up_));
  ntraverse_up_.reference(to_fortran(ntraverse_up_));
  sunpathsfine_up_.reference(to_fortran(sunpathsfine_up_));
  ntraversefine_up_.reference(to_fortran(ntraversefine_up_));
  xfine_dn_.reference(to_fortran(xfine_dn_));
  wfine_dn_.reference(to_fortran(wfine_dn_));
  sunpaths_dn_.reference(to_fortran(sunpaths_dn_));
  ntraverse_dn_.reference(to_fortran(ntraverse_dn_));
  sunpathsfine_dn_.reference(to_fortran(sunpathsfine_dn_));
  ntraversefine_dn_.reference(to_fortran(ntraversefine_dn_));
  xfine_up_p_.reference(to_fortran(xfine_up_p_));
  wfine_up_p_.reference(to_fortran(wfine_up_p_));
  sunpaths_up_p_.reference(to_fortran(sunpaths_up_p_));
  ntraverse_up_p_.reference(to_fortran(ntraverse_up_p_));
  sunpathsfine_up_p_.reference(to_fortran(sunpathsfine_up_p_));
  ntraversefine_up_p_.reference(to_fortran(ntraversefine_up_p_));
  xfine_dn_p_.reference(to_fortran(xfine_dn_p_));
  wfine_dn_p_.reference(to_fortran(wfine_dn_p_));
  sunpaths_dn_p_.reference(to_fortran(sunpaths_dn_p_));
  ntraverse_dn_p_.reference(to_fortran(ntraverse_dn_p_));
  sunpathsfine_dn_p_.reference(to_fortran(sunpathsfine_dn_p_));
  ntraversefine_dn_p_.reference(to_fortran(ntraversefine_dn_p_));
}

////////
// Serialization functions of Fo_Scalarss_Spherfuncs
template<class Archive>
void Fo_Scalarss_Spherfuncs::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Fo_Scalarss_Spherfuncs);
    ar 
       & FP_NVP_(maxmoments)
       & FP_NVP_(maxgeoms)
       & FP_NVP_(nmoments)
       & FP_NVP_(ngeoms)
       & FP_NVP_(starter)
       & FP_NVP_(do_spherfunc)
       & FP_NVP_(cosscat)
       & FP_NVP_(df1)
       & FP_NVP_(df2)
       & FP_NVP_(ss_pleg);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Fo_Scalarss_Spherfuncs::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void Fo_Scalarss_Spherfuncs::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  ss_pleg_.reference(to_fortran(ss_pleg_));
}

////////
// Serialization functions of Fo_Thermal_Rtcalcs_I
template<class Archive>
void Fo_Thermal_Rtcalcs_I::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Fo_Thermal_Rtcalcs_I);
    ar 
       & FP_NVP_(maxgeoms)
       & FP_NVP_(maxlayers)
       & FP_NVP_(maxpartials)
       & FP_NVP_(maxfine)
       & FP_NVP_(max_user_levels)
       & FP_NVP_(do_thermset)
       & FP_NVP_(do_deltam_scaling)
       & FP_NVP_(do_partials)
       & FP_NVP_(do_planpar)
       & FP_NVP_(do_enhanced_ps)
       & FP_NVP_(do_sources_dn)
       & FP_NVP_(do_sources_dn_p)
       & FP_NVP_(ngeoms)
       & FP_NVP_(nlayers)
       & FP_NVP_(nfinedivs)
       & FP_NVP_(n_user_levels)
       & FP_NVP_(user_levels)
       & FP_NVP_(npartials)
       & FP_NVP_(nfinedivs_p)
       & FP_NVP_(partial_outindex)
       & FP_NVP_(partial_outflag)
       & FP_NVP_(partial_layeridx)
       & FP_NVP_(bb_input)
       & FP_NVP_(extinction)
       & FP_NVP_(deltaus)
       & FP_NVP_(omega)
       & FP_NVP_(truncfac)
       & FP_NVP_(mu1)
       & FP_NVP_(losw_paths)
       & FP_NVP_(losp_paths)
       & FP_NVP_(xfine)
       & FP_NVP_(wfine)
       & FP_NVP_(hfine)
       & FP_NVP_(xfine_p)
       & FP_NVP_(wfine_p)
       & FP_NVP_(hfine_p)
       & FP_NVP_(intensity_dta_dn)
       & FP_NVP_(cumsource_dn)
       & FP_NVP_(tcom1)
       & FP_NVP_(do_sources_up)
       & FP_NVP_(do_sources_up_p)
       & FP_NVP_(surfbb)
       & FP_NVP_(user_emissivity)
       & FP_NVP_(intensity_dta_up)
       & FP_NVP_(intensity_dts)
       & FP_NVP_(cumsource_up)
       & FP_NVP_(lostrans_up)
       & FP_NVP_(lostrans_up_p)
       & FP_NVP_(do_upwelling)
       & FP_NVP_(do_dnwelling)
       & FP_NVP_(xfine_up)
       & FP_NVP_(wfine_up)
       & FP_NVP_(hfine_up)
       & FP_NVP_(xfine_dn)
       & FP_NVP_(wfine_dn)
       & FP_NVP_(hfine_dn)
       & FP_NVP_(xfine_up_p)
       & FP_NVP_(wfine_up_p)
       & FP_NVP_(hfine_up_p)
       & FP_NVP_(xfine_dn_p)
       & FP_NVP_(wfine_dn_p)
       & FP_NVP_(hfine_dn_p);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Fo_Thermal_Rtcalcs_I::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void Fo_Thermal_Rtcalcs_I::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  do_sources_dn_.reference(to_fortran(do_sources_dn_));
  do_sources_dn_p_.reference(to_fortran(do_sources_dn_p_));
  nfinedivs_.reference(to_fortran(nfinedivs_));
  nfinedivs_p_.reference(to_fortran(nfinedivs_p_));
  losw_paths_.reference(to_fortran(losw_paths_));
  losp_paths_.reference(to_fortran(losp_paths_));
  xfine_.reference(to_fortran(xfine_));
  wfine_.reference(to_fortran(wfine_));
  hfine_.reference(to_fortran(hfine_));
  xfine_p_.reference(to_fortran(xfine_p_));
  wfine_p_.reference(to_fortran(wfine_p_));
  hfine_p_.reference(to_fortran(hfine_p_));
  intensity_dta_dn_.reference(to_fortran(intensity_dta_dn_));
  cumsource_dn_.reference(to_fortran(cumsource_dn_));
  tcom1_.reference(to_fortran(tcom1_));
  do_sources_up_.reference(to_fortran(do_sources_up_));
  do_sources_up_p_.reference(to_fortran(do_sources_up_p_));
  intensity_dta_up_.reference(to_fortran(intensity_dta_up_));
  intensity_dts_.reference(to_fortran(intensity_dts_));
  cumsource_up_.reference(to_fortran(cumsource_up_));
  lostrans_up_.reference(to_fortran(lostrans_up_));
  lostrans_up_p_.reference(to_fortran(lostrans_up_p_));
  xfine_up_.reference(to_fortran(xfine_up_));
  wfine_up_.reference(to_fortran(wfine_up_));
  hfine_up_.reference(to_fortran(hfine_up_));
  xfine_dn_.reference(to_fortran(xfine_dn_));
  wfine_dn_.reference(to_fortran(wfine_dn_));
  hfine_dn_.reference(to_fortran(hfine_dn_));
  xfine_up_p_.reference(to_fortran(xfine_up_p_));
  wfine_up_p_.reference(to_fortran(wfine_up_p_));
  hfine_up_p_.reference(to_fortran(hfine_up_p_));
  xfine_dn_p_.reference(to_fortran(xfine_dn_p_));
  wfine_dn_p_.reference(to_fortran(wfine_dn_p_));
  hfine_dn_p_.reference(to_fortran(hfine_dn_p_));
}

////////
// Serialization functions of Fo_Thermal_Rtcalcs_Ilps
template<class Archive>
void Fo_Thermal_Rtcalcs_Ilps::serialize(Archive& ar, const unsigned int version)
{
  FP_GENERIC_BASE(Fo_Thermal_Rtcalcs_Ilps);
    ar 
       & FP_NVP_(maxgeoms)
       & FP_NVP_(maxlayers)
       & FP_NVP_(maxpartials)
       & FP_NVP_(maxfine)
       & FP_NVP_(max_user_levels)
       & FP_NVP_(max_atmoswfs)
       & FP_NVP_(do_thermset)
       & FP_NVP_(do_deltam_scaling)
       & FP_NVP_(do_partials)
       & FP_NVP_(do_planpar)
       & FP_NVP_(do_enhanced_ps)
       & FP_NVP_(do_sources_dn)
       & FP_NVP_(do_sources_dn_p)
       & FP_NVP_(do_profilewfs)
       & FP_NVP_(lvaryflags)
       & FP_NVP_(lvarynums)
       & FP_NVP_(ngeoms)
       & FP_NVP_(nlayers)
       & FP_NVP_(nfinedivs)
       & FP_NVP_(n_user_levels)
       & FP_NVP_(user_levels)
       & FP_NVP_(npartials)
       & FP_NVP_(nfinedivs_p)
       & FP_NVP_(partial_outindex)
       & FP_NVP_(partial_outflag)
       & FP_NVP_(partial_layeridx)
       & FP_NVP_(bb_input)
       & FP_NVP_(extinction)
       & FP_NVP_(deltaus)
       & FP_NVP_(omega)
       & FP_NVP_(truncfac)
       & FP_NVP_(l_extinction)
       & FP_NVP_(l_deltaus)
       & FP_NVP_(l_omega)
       & FP_NVP_(l_truncfac)
       & FP_NVP_(mu1)
       & FP_NVP_(losw_paths)
       & FP_NVP_(losp_paths)
       & FP_NVP_(xfine)
       & FP_NVP_(wfine)
       & FP_NVP_(hfine)
       & FP_NVP_(xfine_p)
       & FP_NVP_(wfine_p)
       & FP_NVP_(hfine_p)
       & FP_NVP_(intensity_dta_dn)
       & FP_NVP_(lp_jacobians_dta_dn)
       & FP_NVP_(tcom1)
       & FP_NVP_(l_tcom1)
       & FP_NVP_(max_surfacewfs)
       & FP_NVP_(do_sources_up)
       & FP_NVP_(do_sources_up_p)
       & FP_NVP_(do_surfacewfs)
       & FP_NVP_(n_surfacewfs)
       & FP_NVP_(surfbb)
       & FP_NVP_(user_emissivity)
       & FP_NVP_(ls_user_emissivity)
       & FP_NVP_(intensity_dta_up)
       & FP_NVP_(intensity_dts)
       & FP_NVP_(lp_jacobians_dta_up)
       & FP_NVP_(lp_jacobians_dts_up)
       & FP_NVP_(ls_jacobians_dts)
       & FP_NVP_(lostrans_up)
       & FP_NVP_(lostrans_up_p)
       & FP_NVP_(l_lostrans_up)
       & FP_NVP_(l_lostrans_up_p)
       & FP_NVP_(do_upwelling)
       & FP_NVP_(do_dnwelling)
       & FP_NVP_(xfine_up)
       & FP_NVP_(wfine_up)
       & FP_NVP_(hfine_up)
       & FP_NVP_(xfine_dn)
       & FP_NVP_(wfine_dn)
       & FP_NVP_(hfine_dn)
       & FP_NVP_(xfine_up_p)
       & FP_NVP_(wfine_up_p)
       & FP_NVP_(hfine_up_p)
       & FP_NVP_(xfine_dn_p)
       & FP_NVP_(wfine_dn_p)
       & FP_NVP_(hfine_dn_p);

  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Fo_Thermal_Rtcalcs_Ilps::save(Archive & UNUSED(a), const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void Fo_Thermal_Rtcalcs_Ilps::load(Archive & UNUSED(ar), const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  do_sources_dn_.reference(to_fortran(do_sources_dn_));
  do_sources_dn_p_.reference(to_fortran(do_sources_dn_p_));
  nfinedivs_.reference(to_fortran(nfinedivs_));
  nfinedivs_p_.reference(to_fortran(nfinedivs_p_));
  l_extinction_.reference(to_fortran(l_extinction_));
  l_deltaus_.reference(to_fortran(l_deltaus_));
  l_omega_.reference(to_fortran(l_omega_));
  l_truncfac_.reference(to_fortran(l_truncfac_));
  losw_paths_.reference(to_fortran(losw_paths_));
  losp_paths_.reference(to_fortran(losp_paths_));
  xfine_.reference(to_fortran(xfine_));
  wfine_.reference(to_fortran(wfine_));
  hfine_.reference(to_fortran(hfine_));
  xfine_p_.reference(to_fortran(xfine_p_));
  wfine_p_.reference(to_fortran(wfine_p_));
  hfine_p_.reference(to_fortran(hfine_p_));
  intensity_dta_dn_.reference(to_fortran(intensity_dta_dn_));
  lp_jacobians_dta_dn_.reference(to_fortran(lp_jacobians_dta_dn_));
  tcom1_.reference(to_fortran(tcom1_));
  l_tcom1_.reference(to_fortran(l_tcom1_));
  do_sources_up_.reference(to_fortran(do_sources_up_));
  do_sources_up_p_.reference(to_fortran(do_sources_up_p_));
  ls_user_emissivity_.reference(to_fortran(ls_user_emissivity_));
  intensity_dta_up_.reference(to_fortran(intensity_dta_up_));
  intensity_dts_.reference(to_fortran(intensity_dts_));
  lp_jacobians_dta_up_.reference(to_fortran(lp_jacobians_dta_up_));
  lp_jacobians_dts_up_.reference(to_fortran(lp_jacobians_dts_up_));
  ls_jacobians_dts_.reference(to_fortran(ls_jacobians_dts_));
  lostrans_up_.reference(to_fortran(lostrans_up_));
  lostrans_up_p_.reference(to_fortran(lostrans_up_p_));
  l_lostrans_up_.reference(to_fortran(l_lostrans_up_));
  l_lostrans_up_p_.reference(to_fortran(l_lostrans_up_p_));
  xfine_up_.reference(to_fortran(xfine_up_));
  wfine_up_.reference(to_fortran(wfine_up_));
  hfine_up_.reference(to_fortran(hfine_up_));
  xfine_dn_.reference(to_fortran(xfine_dn_));
  wfine_dn_.reference(to_fortran(wfine_dn_));
  hfine_dn_.reference(to_fortran(hfine_dn_));
  xfine_up_p_.reference(to_fortran(xfine_up_p_));
  wfine_up_p_.reference(to_fortran(wfine_up_p_));
  hfine_up_p_.reference(to_fortran(hfine_up_p_));
  xfine_dn_p_.reference(to_fortran(xfine_dn_p_));
  wfine_dn_p_.reference(to_fortran(wfine_dn_p_));
  hfine_dn_p_.reference(to_fortran(hfine_dn_p_));
}

FP_IMPLEMENT(Fo_Dtwpgeometry_Master);
FP_IMPLEMENT(Fo_Sswpgeometry_Master);
FP_IMPLEMENT(Fo_Scalarss_Rtcalcs_I);
FP_IMPLEMENT(Fo_Scalarss_Rtcalcs_Ilps);
FP_IMPLEMENT(Fo_Scalarss_Spherfuncs);
FP_IMPLEMENT(Fo_Thermal_Rtcalcs_I);
FP_IMPLEMENT(Fo_Thermal_Rtcalcs_Ilps);

#endif