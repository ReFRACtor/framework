#include "first_order_interface.h"
#include "fp_serialize_support.h"
#include "linear_algebra.h"

// This was written by hand. Might be nice to include with automatic
// code generation.

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Fo_Ssgeometry_Master::serialize(Archive& ar,
					     const unsigned int version)
{
  FP_GENERIC_BASE(Fo_Ssgeometry_Master);
  ar & FP_NVP_(maxgeoms)
    & FP_NVP_(maxszas)
    & FP_NVP_(maxvzas)
    & FP_NVP_(maxazms)
    & FP_NVP_(maxlayers)
    & FP_NVP_(maxfine)
    & FP_NVP_(do_obsgeom)
    & FP_NVP_(do_chapman)
    & FP_NVP_(do_planpar)
    & FP_NVP_(do_enhanced_ps)
    & FP_NVP_(ngeoms)
    & FP_NVP_(nszas)
    & FP_NVP_(nvzas)
    & FP_NVP_(nazms)
    & FP_NVP_(nlayers)
    & FP_NVP_(nfine)
    & FP_NVP_(dtr)
    & FP_NVP_(pie)
    & FP_NVP_(vsign)
    & FP_NVP_(eradius)
    & FP_NVP_(heights)
    & FP_NVP_(obsgeom_boa)
    & FP_NVP_(alpha_boa)
    & FP_NVP_(theta_boa)
    & FP_NVP_(phi_boa)
    & FP_NVP_(donadir)
    & FP_NVP_(docrit)
    & FP_NVP_(acrit)
    & FP_NVP_(extinc)
    & FP_NVP_(raycon)
    & FP_NVP_(radii)
    & FP_NVP_(alpha)
    & FP_NVP_(cota)
    & FP_NVP_(nfinedivs)
    & FP_NVP_(xfine)
    & FP_NVP_(wfine)
    & FP_NVP_(csqfine)
    & FP_NVP_(cotfine)
    & FP_NVP_(alphafine)
    & FP_NVP_(radiifine)
    & FP_NVP_(ncrit)
    & FP_NVP_(radcrit)
    & FP_NVP_(cotcrit)
    & FP_NVP_(mu0)
    & FP_NVP_(mu1)
    & FP_NVP_(cosscat)
    & FP_NVP_(chapfacs)
    & FP_NVP_(sunpaths)
    & FP_NVP_(ntraverse)
    & FP_NVP_(sunpathsfine)
    & FP_NVP_(ntraversefine)
    & FP_NVP_(fail)
    & FP_NVP_(message)
    & FP_NVP_(trace);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Fo_Ssgeometry_Master::save(Archive & UNUSED(a),
		    const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}

template<class Archive>
void Fo_Ssgeometry_Master::load(Archive & UNUSED(ar),
			    const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  obsgeom_boa_.reference(to_fortran(obsgeom_boa_));
  alpha_.reference(to_fortran(alpha_));
  cota_.reference(to_fortran(cota_));
  nfinedivs_.reference(to_fortran(nfinedivs_));
  xfine_.reference(to_fortran(xfine_));
  wfine_.reference(to_fortran(wfine_));
  csqfine_.reference(to_fortran(csqfine_));
  cotfine_.reference(to_fortran(cotfine_));
  alphafine_.reference(to_fortran(alphafine_));
  radiifine_.reference(to_fortran(radiifine_));
  chapfacs_.reference(to_fortran(chapfacs_));
  sunpaths_.reference(to_fortran(sunpaths_));
  ntraverse_.reference(to_fortran(ntraverse_));
  sunpathsfine_.reference(to_fortran(sunpathsfine_));
  ntraversefine_.reference(to_fortran(ntraversefine_));
}

template<class Archive>
void Fo_Scalarss_Spherfuncs::serialize(Archive& ar,
					     const unsigned int version)
{
  FP_GENERIC_BASE(Fo_Scalarss_Spherfuncs);
  ar & FP_NVP_(starter)
    & FP_NVP_(maxmoms)
    & FP_NVP_(maxgeoms)
    & FP_NVP_(nmoms)
    & FP_NVP_(ngeoms)
    & FP_NVP_(df1)
    & FP_NVP_(df2)
    & FP_NVP_(cosscat)
    & FP_NVP_(ss_pleg);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Fo_Scalarss_Spherfuncs::save(Archive & UNUSED(a),
		    const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}
template<class Archive>
void Fo_Scalarss_Spherfuncs::load(Archive & UNUSED(ar),
			    const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  ss_pleg_.reference(to_fortran(ss_pleg_));
}

template<class Archive>
void Fo_Scalarss_Rtcalcs_Ilps::serialize(Archive& ar,
					     const unsigned int version)
{
  FP_GENERIC_BASE(Fo_Scalarss_Rtcalcs_Ilps);
  ar & FP_NVP_(maxgeoms)
    & FP_NVP_(maxlayers)
    & FP_NVP_(maxfine)
    & FP_NVP_(max_atmoswfs)
    & FP_NVP_(max_surfacewfs)
    & FP_NVP_(do_planpar)
    & FP_NVP_(do_regular_ps)
    & FP_NVP_(do_enhanced_ps)
    & FP_NVP_(donadir)
    & FP_NVP_(do_sleave)
    & FP_NVP_(do_profilewfs)
    & FP_NVP_(do_reflecwfs)
    & FP_NVP_(do_sleavewfs)
    & FP_NVP_(lvaryflags)
    & FP_NVP_(lvarynums)
    & FP_NVP_(n_reflecwfs)
    & FP_NVP_(n_sleavewfs)
    & FP_NVP_(ngeoms)
    & FP_NVP_(nlayers)
    & FP_NVP_(nfinedivs)
    & FP_NVP_(aclevel)
    & FP_NVP_(reflec)
    & FP_NVP_(slterm)
    & FP_NVP_(extinction)
    & FP_NVP_(deltaus)
    & FP_NVP_(exactscat_up)
    & FP_NVP_(flux)
    & FP_NVP_(ls_reflec)
    & FP_NVP_(lssl_slterm)
    & FP_NVP_(l_extinction)
    & FP_NVP_(l_deltaus)
    & FP_NVP_(l_exactscat_up)
    & FP_NVP_(mu0)
    & FP_NVP_(mu1)
    & FP_NVP_(ncrit)
    & FP_NVP_(xfine)
    & FP_NVP_(wfine)
    & FP_NVP_(csqfine)
    & FP_NVP_(cotfine)
    & FP_NVP_(raycon)
    & FP_NVP_(cota)
    & FP_NVP_(sunpaths)
    & FP_NVP_(ntraverse)
    & FP_NVP_(sunpaths_fine)
    & FP_NVP_(ntraverse_fine)
    & FP_NVP_(intensity_up)
    & FP_NVP_(intensity_db)
    & FP_NVP_(lp_jacobians_up)
    & FP_NVP_(lp_jacobians_db)
    & FP_NVP_(ls_jacobians_db);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Fo_Scalarss_Rtcalcs_Ilps::save(Archive & UNUSED(a),
		    const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}
template<class Archive>
void Fo_Scalarss_Rtcalcs_Ilps::load(Archive & UNUSED(ar),
			    const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  nfinedivs_.reference(to_fortran(nfinedivs_));
  exactscat_up_.reference(to_fortran(exactscat_up_));
  ls_reflec_.reference(to_fortran(ls_reflec_));
  lssl_slterm_.reference(to_fortran(lssl_slterm_));
  l_extinction_.reference(to_fortran(l_extinction_));
  l_deltaus_.reference(to_fortran(l_deltaus_));
  l_exactscat_up_.reference(to_fortran(l_exactscat_up_));
  xfine_.reference(to_fortran(xfine_));
  wfine_.reference(to_fortran(wfine_));
  csqfine_.reference(to_fortran(csqfine_));
  cotfine_.reference(to_fortran(cotfine_));
  cota_.reference(to_fortran(cota_));
  sunpaths_.reference(to_fortran(sunpaths_));
  ntraverse_.reference(to_fortran(ntraverse_));
  sunpaths_fine_.reference(to_fortran(sunpaths_fine_));
  ntraverse_fine_.reference(to_fortran(ntraverse_fine_));
  lp_jacobians_up_.reference(to_fortran(lp_jacobians_up_));
  lp_jacobians_db_.reference(to_fortran(lp_jacobians_db_));
  ls_jacobians_db_.reference(to_fortran(ls_jacobians_db_));
}

FP_IMPLEMENT(Fo_Ssgeometry_Master);
FP_IMPLEMENT(Fo_Scalarss_Spherfuncs);
FP_IMPLEMENT(Fo_Scalarss_Rtcalcs_Ilps);
#endif
