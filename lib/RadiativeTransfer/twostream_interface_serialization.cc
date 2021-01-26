#include "twostream_interface.h"
#include "fp_serialize_support.h"
#include "linear_algebra.h"

// This was written by hand. Might be nice to include with automatic
// code generation.

using namespace FullPhysics;

#ifdef FP_HAVE_BOOST_SERIALIZATION
template<class Archive>
void Twostream_Ls_Brdf_Supplement::serialize(Archive& ar,
					     const unsigned int version)
{
  FP_GENERIC_BASE(Twostream_Ls_Brdf_Supplement);
  ar & FP_NVP_(maxbeams)
    & FP_NVP_(max_user_streams)
    & FP_NVP_(max_user_obsgeoms)
    & FP_NVP_(maxstreams_brdf)
    & FP_NVP_(max_brdf_kernels)
    & FP_NVP_(max_brdf_parameters)
    & FP_NVP_(max_surfacewfs)
    & FP_NVP_(do_solar_sources)
    & FP_NVP_(do_user_obsgeoms)
    & FP_NVP_(lambertian_kernel_flag)
    & FP_NVP_(do_shadow_effect)
    & FP_NVP_(do_surface_emission)
    & FP_NVP_(nbeams)
    & FP_NVP_(n_user_streams)
    & FP_NVP_(n_user_obsgeoms)
    & FP_NVP_(beam_szas)
    & FP_NVP_(user_angles)
    & FP_NVP_(user_obsgeoms)
    & FP_NVP_(stream_value)
    & FP_NVP_(nstreams_brdf)
    & FP_NVP_(n_brdf_kernels)
    & FP_NVP_(which_brdf)
    & FP_NVP_(brdf_factors)
    & FP_NVP_(n_brdf_parameters)
    & FP_NVP_(brdf_parameters)
    & FP_NVP_(do_kernel_factor_wfs)
    & FP_NVP_(do_kernel_params_wfs)
    & FP_NVP_(do_kparams_derivs)
    & FP_NVP_(n_surface_wfs)
    & FP_NVP_(n_kernel_factor_wfs)
    & FP_NVP_(n_kernel_params_wfs)
    & FP_NVP_(brdf_f_0)
    & FP_NVP_(brdf_f)
    & FP_NVP_(ubrdf_f)
    & FP_NVP_(emissivity)
    & FP_NVP_(ls_brdf_f_0)
    & FP_NVP_(ls_brdf_f)
    & FP_NVP_(ls_ubrdf_f)
    & FP_NVP_(ls_emissivity)
    & FP_NVP_(status_brdfsup)
    & FP_NVP_(message)
    & FP_NVP_(action);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Twostream_Ls_Brdf_Supplement::save(Archive & UNUSED(a),
		    const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}
template<class Archive>
void Twostream_Ls_Brdf_Supplement::load(Archive & UNUSED(ar),
			    const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  user_obsgeoms_.reference(to_fortran(user_obsgeoms_));
  brdf_parameters_.reference(to_fortran(brdf_parameters_));
  do_kernel_params_wfs_.reference(to_fortran(do_kernel_params_wfs_));
  brdf_f_0_.reference(to_fortran(brdf_f_0_));
  ubrdf_f_.reference(to_fortran(ubrdf_f_));
  ls_brdf_f_0_.reference(to_fortran(ls_brdf_f_0_));
  ls_brdf_f_.reference(to_fortran(ls_brdf_f_));
  ls_ubrdf_f_.reference(to_fortran(ls_ubrdf_f_));
}

template<class Archive>
void Twostream_Lps_Master::serialize(Archive& ar,
					     const unsigned int version)
{
  FP_GENERIC_BASE(Twostream_Lps_Master);
  ar & FP_NVP_(maxlayers)
    & FP_NVP_(maxtotal)
    & FP_NVP_(maxmessages)
    & FP_NVP_(maxbeams)
    & FP_NVP_(max_geometries)
    & FP_NVP_(max_user_streams)
    & FP_NVP_(max_user_relazms)
    & FP_NVP_(max_user_obsgeoms)
    & FP_NVP_(max_atmoswfs)
    & FP_NVP_(max_surfacewfs)
    & FP_NVP_(max_sleavewfs)
    & FP_NVP_(do_upwelling)
    & FP_NVP_(do_dnwelling)
    & FP_NVP_(do_plane_parallel)
    & FP_NVP_(do_2s_levelout)
    & FP_NVP_(do_mvout_only)
    & FP_NVP_(do_additional_mvout)
    & FP_NVP_(do_solar_sources)
    & FP_NVP_(do_thermal_emission)
    & FP_NVP_(do_surface_emission)
    & FP_NVP_(do_d2s_scaling)
    & FP_NVP_(do_brdf_surface)
    & FP_NVP_(do_user_obsgeoms)
    & FP_NVP_(do_surface_leaving)
    & FP_NVP_(do_sl_isotropic)
    & FP_NVP_(do_pentadiag_inverse)
    & FP_NVP_(bvpindex)
    & FP_NVP_(bvpscalefactor)
    & FP_NVP_(taylor_order)
    & FP_NVP_(taylor_small)
    & FP_NVP_(tcutoff)
    & FP_NVP_(nlayers)
    & FP_NVP_(ntotal)
    & FP_NVP_(stream_value)
    & FP_NVP_(n_user_obsgeoms)
    & FP_NVP_(user_obsgeoms)
    & FP_NVP_(n_user_streams)
    & FP_NVP_(user_angles)
    & FP_NVP_(n_user_relazms)
    & FP_NVP_(user_relazms)
    & FP_NVP_(flux_factor)
    & FP_NVP_(nbeams)
    & FP_NVP_(beam_szas)
    & FP_NVP_(earth_radius)
    & FP_NVP_(height_grid)
    & FP_NVP_(deltau_input)
    & FP_NVP_(omega_input)
    & FP_NVP_(asymm_input)
    & FP_NVP_(d2s_scaling)
    & FP_NVP_(thermal_bb_input)
    & FP_NVP_(lambertian_albedo)
    & FP_NVP_(brdf_f_0)
    & FP_NVP_(brdf_f)
    & FP_NVP_(ubrdf_f)
    & FP_NVP_(emissivity)
    & FP_NVP_(surfbb)
    & FP_NVP_(slterm_isotropic)
    & FP_NVP_(slterm_f_0)
    & FP_NVP_(do_profile_wfs)
    & FP_NVP_(do_surface_wfs)
    & FP_NVP_(do_sleave_wfs)
    & FP_NVP_(layer_vary_flag)
    & FP_NVP_(layer_vary_number)
    & FP_NVP_(n_surface_wfs)
    & FP_NVP_(n_sleave_wfs)
    & FP_NVP_(lssl_slterm_isotropic)
    & FP_NVP_(lssl_slterm_f_0)
    & FP_NVP_(l_deltau_input)
    & FP_NVP_(l_omega_input)
    & FP_NVP_(l_asymm_input)
    & FP_NVP_(l_d2s_scaling)
    & FP_NVP_(ls_brdf_f_0)
    & FP_NVP_(ls_brdf_f)
    & FP_NVP_(ls_ubrdf_f)
    & FP_NVP_(ls_emissivity)
    & FP_NVP_(intensity_toa)
    & FP_NVP_(profilewf_toa)
    & FP_NVP_(surfacewf_toa)
    & FP_NVP_(intensity_boa)
    & FP_NVP_(profilewf_boa)
    & FP_NVP_(surfacewf_boa)
    & FP_NVP_(radlevel_up)
    & FP_NVP_(radlevel_dn)
    & FP_NVP_(n_geometries)
    & FP_NVP_(profjaclevel_up)
    & FP_NVP_(profjaclevel_dn)
    & FP_NVP_(surfjaclevel_up)
    & FP_NVP_(surfjaclevel_dn)
    & FP_NVP_(fluxes_toa)
    & FP_NVP_(profjacfluxes_toa)
    & FP_NVP_(surfjacfluxes_toa)
    & FP_NVP_(fluxes_boa)
    & FP_NVP_(profjacfluxes_boa)
    & FP_NVP_(surfjacfluxes_boa)
    & FP_NVP_(status_inputcheck)
    & FP_NVP_(c_nmessages)
    & FP_NVP_(c_messages)
    & FP_NVP_(c_actions)
    & FP_NVP_(status_execution)
    & FP_NVP_(e_message)
    & FP_NVP_(e_trace_1)
    & FP_NVP_(e_trace_2);
  boost::serialization::split_member(ar, *this, version);
}

template<class Archive>
void Twostream_Lps_Master::save(Archive & UNUSED(a),
		    const unsigned int UNUSED(version)) const
{
  // Nothing more to do
}
template<class Archive>
void Twostream_Lps_Master::load(Archive & UNUSED(ar),
			    const unsigned int UNUSED(version))
{
  // The arrays come back from serialization in C order. We need them
  // in fortran order, so change them to that.
  user_obsgeoms_.reference(to_fortran(user_obsgeoms_));
  brdf_f_0_.reference(to_fortran(brdf_f_0_));
  ubrdf_f_.reference(to_fortran(ubrdf_f_));
  slterm_f_0_.reference(to_fortran(slterm_f_0_));
  lssl_slterm_isotropic_.reference(to_fortran(lssl_slterm_isotropic_));
  lssl_slterm_f_0_.reference(to_fortran(lssl_slterm_f_0_));
  l_deltau_input_.reference(to_fortran(l_deltau_input_));
  l_omega_input_.reference(to_fortran(l_omega_input_));
  l_asymm_input_.reference(to_fortran(l_asymm_input_));
  l_d2s_scaling_.reference(to_fortran(l_d2s_scaling_));
  ls_brdf_f_0_.reference(to_fortran(ls_brdf_f_0_));
  ls_brdf_f_.reference(to_fortran(ls_brdf_f_));
  ls_ubrdf_f_.reference(to_fortran(ls_ubrdf_f_));
  profilewf_toa_.reference(to_fortran(profilewf_toa_));
  surfacewf_toa_.reference(to_fortran(surfacewf_toa_));
  profilewf_boa_.reference(to_fortran(profilewf_boa_));
  surfacewf_boa_.reference(to_fortran(surfacewf_boa_));
  radlevel_up_.reference(to_fortran(radlevel_up_));
  radlevel_dn_.reference(to_fortran(radlevel_dn_));
  profjaclevel_up_.reference(to_fortran(profjaclevel_up_));
  profjaclevel_dn_.reference(to_fortran(profjaclevel_dn_));
  surfjaclevel_up_.reference(to_fortran(surfjaclevel_up_));
  surfjaclevel_dn_.reference(to_fortran(surfjaclevel_dn_));
  fluxes_toa_.reference(to_fortran(fluxes_toa_));
  profjacfluxes_toa_.reference(to_fortran(profjacfluxes_toa_));
  surfjacfluxes_toa_.reference(to_fortran(surfjacfluxes_toa_));
  fluxes_boa_.reference(to_fortran(fluxes_boa_));
  profjacfluxes_boa_.reference(to_fortran(profjacfluxes_boa_));
  surfjacfluxes_boa_.reference(to_fortran(surfjacfluxes_boa_));
  c_messages_.reference(to_fortran(c_messages_));
  c_actions_.reference(to_fortran(c_actions_));
}

FP_IMPLEMENT(Twostream_Ls_Brdf_Supplement);
FP_IMPLEMENT(Twostream_Lps_Master);
#endif
